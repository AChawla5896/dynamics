BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/sources.mk

.SECONDARY: $(simp_tests_OBJS)
TEST=$(BLD_DIR)/simp/tests/Test

all: $(simp_tests_EXES)

run: $(TEST)
	$(TEST) $(SRC_DIR)/simp/tests/ 
	cd $(SRC_DIR)/simp/tests; $(MAKE) clean-outputs

quiet: $(TEST)
	$(TEST) $(SRC_DIR)/simp/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in simp/tests/log" > count
	@cat count
	@cd $(SRC_DIR)/simp/tests; $(MAKE) clean-outputs

clean:
	rm -f $(simp_tests_EXES)
	rm -f $(simp_tests_OBJS) 
	rm -f $(simp_tests_OBJS:.o=.d) 
	rm -f log count
ifeq ($(BLD_DIR),$(SRC_DIR))
	cd analysis; $(MAKE) clean
	cd boundary; $(MAKE) clean
	cd interaction; $(MAKE) clean
	cd species; $(MAKE) clean
else
	cd $(SRC_DIR)/simp/tests; $(MAKE) clean-outputs
endif

clean-outputs:
	@cd species; $(MAKE) clean-outputs
	@cd analysis; $(MAKE) clean-outputs

   

-include $(simp_tests_OBJS:.o=.d)
