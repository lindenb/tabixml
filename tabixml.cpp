#include <libxml/xmlreader.h>
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xmlIO.h>
#include "bgzf.h"
#include <ctype.h>
#include <assert.h>
#include <sys/stat.h>
#include "khash.h"
#include "ksort.h"
#include "kstring.h"
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

using namespace std;

#define TAD_MIN_CHUNK_GAP 32768
// 1<<14 is the size of minimum bin.
#define TAD_LIDX_SHIFT    14

typedef struct {
	uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(offt, pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} ti_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} ti_lidx_t;

KHASH_MAP_INIT_INT(i, ti_binlist_t)
KHASH_MAP_INIT_STR(s, int)

struct __ti_index_t {
	ti_conf_t conf;
	int32_t n, max;
	khash_t(s) *tname;
	khash_t(i) **index;
	ti_lidx_t *index2;
};

struct __ti_iter_t {
	int from_first; // read from the first record; no random access
	int tid, beg, end, n_off, i, finished;
	uint64_t curr_off;
	kstring_t str;
	const ti_index_t *idx;
	pair64_t *off;
};

typedef struct {
	int tid, beg, end, bin;
} ti_intv_t;

static int _xmlOutputWriteCallbackBGZF(
     void * context,
     const char * buffer,
     int len)
    {
    if(context==0) return  len;
    return bgzf_write(((BGZF*)(context)),(const void*)buffer,len);
    }

static int _xmlOutputCloseCallbackBGZF(void * context)
    {
    return bgzf_close(((BGZF*)(context)));
    }
   
xmlOutputBufferPtr xmlOutputBufferCreateBGZF(
	BGZF* out,
	xmlCharEncodingHandlerPtr encoder
	)
    {
    xmlOutputBufferPtr outbuf=::xmlOutputBufferCreateIO(
	_xmlOutputWriteCallbackBGZF,
	_xmlOutputCloseCallbackBGZF,
	 (void*)out,
	 encoder);
    return outbuf;
    }

xmlNodePtr read_next_node(xmlTextReaderPtr in)
	{
	return 0;
	}
static int get_intv(ti_index_t *idx, xmlNodePtr node, ti_intv_t *intv)
	{
	return -1;
	}
static inline void insert_offset(khash_t(i) *h, int bin, uint64_t beg, uint64_t end)
{
	khint_t k;
	ti_binlist_t *l;
	int ret;
	k = kh_put(i, h, bin, &ret);
	l = &kh_value(h, k);
	if (ret) { // not present
		l->m = 1; l->n = 0;
		l->list = (pair64_t*)calloc(l->m, 16);
	}
	if (l->n == l->m) {
		l->m <<= 1;
		l->list = (pair64_t*)realloc(l->list, l->m * 16);
	}
	l->list[l->n].u = beg; l->list[l->n++].v = end;
}

static inline uint64_t insert_offset2(ti_lidx_t *index2, int _beg, int _end, uint64_t offset)
{
	int i, beg, end;
	beg = _beg >> TAD_LIDX_SHIFT;
	end = (_end - 1) >> TAD_LIDX_SHIFT;
	if (index2->m < end + 1) {
		int old_m = index2->m;
		index2->m = end + 1;
		kroundup32(index2->m);
		index2->offset = (uint64_t*)realloc(index2->offset, index2->m * 8);
		memset(index2->offset + old_m, 0, 8 * (index2->m - old_m));
	}
	if (beg == end) {
		if (index2->offset[beg] == 0) index2->offset[beg] = offset;
	} else {
		for (i = beg; i <= end; ++i)
			if (index2->offset[i] == 0) index2->offset[i] = offset;
	}
	if (index2->n < end + 1) index2->n = end + 1;
	return (uint64_t)beg<<32 | end;
}


static void merge_chunks(ti_index_t *idx)
{
	khash_t(i) *index;
	int i, l, m;
	khint_t k;
	for (i = 0; i < idx->n; ++i) {
		index = idx->index[i];
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			ti_binlist_t *p;
			if (!kh_exist(index, k)) continue;
			p = &kh_value(index, k);
			m = 0;
			for (l = 1; l < p->n; ++l) {
				if (p->list[m].v>>16 == p->list[l].u>>16) p->list[m].v = p->list[l].v;
				else p->list[++m] = p->list[l];
			} // ~for(l)
			p->n = m + 1;
		} // ~for(k)
	} // ~for(i)
}

static void fill_missing(ti_index_t *idx)
{
	int i, j;
	for (i = 0; i < idx->n; ++i) {
		ti_lidx_t *idx2 = &idx->index2[i];
		for (j = 1; j < idx2->n; ++j)
			if (idx2->offset[j] == 0)
				idx2->offset[j] = idx2->offset[j-1];
	}
}


ti_index_t *xml_ti_index_core(xmlTextReaderPtr fp,BGZF *out, const ti_conf_t *conf)
	{
	ti_index_t *idx;
	uint32_t last_bin, save_bin;
	int32_t last_coor, last_tid, save_tid;
	uint64_t save_off, last_off, offset0 = (uint64_t)-1, tmp;
	//kstring_t *str;
	xmlNodePtr node;
	//str = (kstring_t *)calloc(1, sizeof(kstring_t));

	idx = (ti_index_t*)calloc(1, sizeof(ti_index_t));
	idx->conf = *conf;
	idx->n = idx->max = 0;
	idx->tname = kh_init(s);
	idx->index = 0;
	idx->index2 = 0;

	save_bin = save_tid = last_tid = last_bin = 0xffffffffu;
	save_off = last_off = bgzf_tell(out); last_coor = 0xffffffffu;
	
	for(;;)
		{
		int ret = xmlTextReaderRead(fp);
		if(ret<=0) break;
		int nodeType=xmlTextReaderNodeType(reader);
    		switch(nodeType)
		    {
		    case XML_READER_TYPE_ELEMENT:
			{
			break;
			}
		    default:break;
		    }
            		
		}
	
	while ((node =read_next_node(fp))!=0)
		{
		ti_intv_t intv;
		get_intv(idx, node, &intv);
		if(intv.beg<0 || intv.end<0 )
			{
		    	fprintf(stderr,"[ti_index_core] the indexes overlap or are out of bounds\n");
		    	exit(EXIT_FAILURE);
			}
		if (last_tid != intv.tid)
			{
			// change of chromosomes
            if (last_tid>intv.tid )
            	{
                fprintf(stderr,"[ti_index_core] the chromosome blocks not continuous at is the file sorted? [pos %d]\n",intv.beg+1);
                exit(EXIT_FAILURE);
            	}
			last_tid = intv.tid;
			last_bin = 0xffffffffu;
		} else if (last_coor > intv.beg) {
			fprintf(stderr, "[ti_index_core] the file out of order at line %llu\n", (unsigned long long)lineno);
			exit(1);
		}
		tmp = insert_offset2(&idx->index2[intv.tid], intv.beg, intv.end, last_off);
		if (last_off == 0) offset0 = tmp;
		if (intv.bin != last_bin) { // then possibly write the binning index
			if (save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
				insert_offset(idx->index[save_tid], save_bin, save_off, last_off);
			save_off = last_off;
			save_bin = last_bin = intv.bin;
			save_tid = intv.tid;
			if (save_tid < 0) break;
		}
		if (bgzf_tell(out) <= last_off) {
			fprintf(stderr, "[ti_index_core] bug in BGZF: %llx < %llx\n",
					(unsigned long long)bgzf_tell(out), (unsigned long long)last_off);
			exit(1);
		}
		last_off = bgzf_tell(out);
		last_coor = intv.beg;
	}
	if (save_tid >= 0) insert_offset(idx->index[save_tid], save_bin, save_off, bgzf_tell(out));
	merge_chunks(idx);
	fill_missing(idx);
	if (offset0 != (uint64_t)-1 && idx->n && idx->index2[0].offset) {
		int i, beg = offset0>>32, end = offset0&0xffffffffu;
		for (i = beg; i <= end; ++i) idx->index2[0].offset[i] = 0;
	}


	//TODO free node ?
	return idx;
	}



int main(int argc,char** argv)
	{
	LIBXML_TEST_VERSION;
	xmlTextReaderPtr xmlReader=xmlReaderForFd(fileno(stdin),0,"UTF-8",0);
	
	xml_ti_index_core(0,0,0);
	return 0;
	}
