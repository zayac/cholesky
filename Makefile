TARGETS    = single openmp symposdef decompose_block decompose
DEPS       = 
SNETC	   = snetc
SNETCFLAGS = -v1 -g -O3 -lm -q -distrib nodist -threading front
CLEANS     = 
CFLAGS     = -Wall -g -O3 

.PRECIOUS: single symposdef

targets: $(TARGETS)

opt:
	$(CC) $(CFLAGS) -O3 single.c -o single-opt -lm
	$(CC) $(CFLAGS) -O3 symposdef.c -o symposdef-opt -lm

single: single.c
	$(CC) $(CFLAGS) -o $@ $< -lm -lrt
	./$@

openmp: openmp.c
	$(CC) -fopenmp $(CFLAGS) -o $@ $< -lm -lrt
	./$@

symposdef: symposdef.c
	$(CC) $(CFLAGS) -o $@ $< -lm
	./$@

decompose_block: decompose_block.snet deblockboxes.o
	$(SNETC) $(SNETCFLAGS) $^ -lrt -o $@

decompose: decompose.snet deboxes.o
	$(SNETC) $(SNETCFLAGS) $^ -lrt -o $@

deblockboxes.o: deblockboxes.c
	$(CC) $(CFLAGS) -c -I$(SNET_INCLUDES) $<

deboxes.o: deboxes.c
	$(CC) $(CFLAGS) -c -I$(SNET_INCLUDES) $<

tags: single.c symposdef.c deblockboxes.c deboxes.c
	ctags $^

clean:
	$(RM) $(RMFLAGS) -- *.o *.a *.lo *.la *.Plo $(TARGETS) core vgcore.*
	$(RM) $(RMFLAGS) -- $(patsubst %.snet,%.[ch],decompose.snet decompose_block.snet)

