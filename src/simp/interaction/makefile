#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(simp_interaction_OBJS)

clean:
	rm -f $(simp_OBJS) $(simp_OBJS:.o=.d) $(simp_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	rm -f */*.o */*.d
	rm -f lib*.a

#-----------------------------------------------------------------------
# Include dependency files

-include $(simp_OBJS:.o=.d)

