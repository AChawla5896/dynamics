.PHONY: print

BLD_DIR_REL = ../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk


USR_DIR = $(SRC_DIR)/mcMd/user

dynamics_prog = $(USR_DIR)/dynamics
cluster_Class = $(USR_DIR)/clusterInfo
mapping_Class = $(USR_DIR)/mapping

analyze_summary_prog = $(USR_DIR)/analyzeSummary
readSample_Class = $(USR_DIR)/readSample

LOCAL_OBJ= $(dynamics_prog).o $(cluster_Class).o $(mapping_Class).o

print:
	@echo $(INCLUDES) $(CXXFLAGS) $(SRC_DIR) $(user_EXE) $(MCMD_DEFINES)

all: $(LOCAL_OBJ) $(dynamics_prog) $(analyze_summary_prog).o $(analyze_summary_prog)

$(dynamics_prog): $(LIBS)
	$(CXX) -o $(dynamics_prog) $(LOCAL_OBJ) $(LIBS) $(LDFLAGS)

$(analyze_summary_prog): $(LIBS)
	$(CXX) -o $(analyze_summary_prog) $(analyze_summary_prog).o $(LIBS) $(LDFLAGS)

clean:
	rm -f $(mcMd_user_OBJS)
	rm -f $(mcMd_user_OBJS:.o=.d)
	rm -f $(LOCAL_OBJ)
	rm -f $(LOCAL_OBJ:.o=.d)
	rm -f $(dynamics_prog)
	rm -f $(dynamics_prog).o
	rm -f $(dynamics_prog).d
	rm -f $(analyze_summary_prog)
	rm -f $(analyze_summary_prog).o
	rm -f $(analyze_summary_prog).d

clean-deps:
	rm -f $(mcMd_user_OBJS:.o=.d)


-include $(mcMd_user_OBJS:.o=.d)
-include $(LOCAL_OBJ:.o=.d)
