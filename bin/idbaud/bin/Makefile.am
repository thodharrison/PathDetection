AM_LDFLAGS = -fopenmp -pthread
LIBS = $(top_srcdir)/lib/libassembly.a
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

sort_reads_SOURCES = $(top_srcdir)/src/tools/sort_reads.cpp
parallel_rna_blat_SOURCES = $(top_srcdir)/src/tools/parallel_rna_blat.cpp
test_SOURCES = $(top_srcdir)/src/tools/test.cpp
print_graph_SOURCES = $(top_srcdir)/src/tools/print_graph.cpp
validate_rna_SOURCES = $(top_srcdir)/src/tools/validate_rna.cpp
scaffold_SOURCES = $(top_srcdir)/src/tools/scaffold.cpp
shuffle_reads_SOURCES = $(top_srcdir)/src/tools/shuffle_reads.cpp
fa2fq_SOURCES = $(top_srcdir)/src/tools/fa2fq.cpp
validate_contigs_mummer_SOURCES = $(top_srcdir)/src/tools/validate_contigs_mummer.cpp
filter_blat_SOURCES = $(top_srcdir)/src/tools/filter_blat.cpp
sort_psl_SOURCES = $(top_srcdir)/src/tools/sort_psl.cpp
validate_reads_blat_SOURCES = $(top_srcdir)/src/tools/validate_reads_blat.cpp
validate_component_SOURCES = $(top_srcdir)/src/tools/validate_component.cpp
filter_contigs_SOURCES = $(top_srcdir)/src/tools/filter_contigs.cpp
sample_reads_SOURCES = $(top_srcdir)/src/tools/sample_reads.cpp
split_fq_SOURCES = $(top_srcdir)/src/tools/split_fq.cpp
idba_tran_test_SOURCES = $(top_srcdir)/src/tools/idba_tran_test.cpp
split_fa_SOURCES = $(top_srcdir)/src/tools/split_fa.cpp 

idba_hybrid_SOURCES = $(top_srcdir)/src/release/idba_hybrid.cpp
fq2fa_SOURCES = $(top_srcdir)/src/release/fq2fa.cpp
raw_n50_SOURCES = $(top_srcdir)/src/release/raw_n50.cpp
parallel_blat_SOURCES = $(top_srcdir)/src/release/parallel_blat.cpp
split_scaffold_SOURCES = $(top_srcdir)/src/release/split_scaffold.cpp
validate_contigs_blat_SOURCES = $(top_srcdir)/src/release/validate_contigs_blat.cpp
filterfa_SOURCES = $(top_srcdir)/src/release/filterfa.cpp
sim_reads_SOURCES = $(top_srcdir)/src/release/sim_reads.cpp
idba_ud_SOURCES = $(top_srcdir)/src/release/idba_ud.cpp
sim_reads_tran_SOURCES = $(top_srcdir)/src/release/sim_reads_tran.cpp
idba_tran_SOURCES = $(top_srcdir)/src/release/idba_tran.cpp
idba_SOURCES = $(top_srcdir)/src/release/idba.cpp 

bin_PROGRAMS = \
	idba_hybrid 

noinst_PROGRAMS = \
	sort_reads \
	parallel_rna_blat \
	test \
	print_graph \
	validate_rna \
	scaffold \
	shuffle_reads \
	fa2fq \
	validate_contigs_mummer \
	filter_blat \
	sort_psl \
	validate_reads_blat \
	validate_component \
	filter_contigs \
	sample_reads \
	split_fq \
	idba_tran_test \
	split_fa \
	fq2fa \
	raw_n50 \
	parallel_blat \
	split_scaffold \
	validate_contigs_blat \
	filterfa \
	sim_reads \
	idba_ud \
	sim_reads_tran \
	idba_tran \
	idba 

