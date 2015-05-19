CC=g++
CFLAGS=-O3 -Wall
htslib.version=1.2.1

LIBS=
CFLAGS+=`xml2-config --cflags`
LIBS+=`xml2-config --libs`
CFLAGS+= -Ihtslib-${htslib.version}/htslib
LIBS+= -Lhtslib-${htslib.version}/ 

tabixml: tabixml.cpp htslib-${htslib.version}/libhts.a
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

htslib-${htslib.version}/libhts.a:
	rm -rf htslib-${htslib.version}
	wget -O ${htslib.version}.tar.gz "https://github.com/samtools/htslib/archive/${htslib.version}.tar.gz"
	tar xvfz ${htslib.version}.tar.gz
	rm ${htslib.version}.tar.gz
	(cd  htslib-${htslib.version} && make libhts.a)

clean:
	rm -f tabixml htslib-${htslib.version} *.o 
