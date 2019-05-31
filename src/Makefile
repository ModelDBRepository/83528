CC = g++ -Wall -O3
CPP = g++ -Wall -O3
LDLIBS = -lm # -lefence
LEX  = flex -I
YACC = byacc -d -v -t

USE_BTTV   = 0
USE_GAZEBO = 0
USE_MIRO   = 0
USE_MIRROR = 0

TARGET     = ../cococo

# measure time: compile: -pg; link: -ldl -L/usr/lib/libp; -lm behind -ll -ly !
# ... compile ... run .... gprof programname | less

PARSER_OBJECTS      = parser/y.tab.o parser/lex.yy.o
KERNEL_OBJECTS      = kernel/coco.o kernel/utils.o kernel/vehicle.o kernel/relax.o
FUN_BASIC_OBJECTS   = fun_basic/local.o fun_weight_list/weight_list.o fun_weight_list/weight_sigpi.o fun_basic/observe.o fun_basic/single.o fun_basic/total.o
FUN_SPECIAL_OBJECTS = fun_special/reinforcement.o
FUN_DATA_OBJECTS    = fun_data/data.made.o fun_data/data.images.o
                      # fun_data/data.digits.o fun_data/data.motor.o

ifeq ($(USE_BTTV),1)
BTTVOBJECTS        = fun_data/data.bttv.o
endif

ifeq ($(USE_GAZEBO),1)
FUN_GAZEBO_OBJECTS   = fun_gazebo/gazebo.io.o
GAZEBOFLAGS          = -g3 `pkg-config --cflags gazebo`
GAZEBOLIBS           = `pkg-config --libs gazebo`
# gcc -Wall -g3 `pkg-config --cflags gazebo` -c -o simple.o simple.c
# gcc simple.o -o simple `pkg-config --libs gazebo`
endif

ifeq ($(USE_MIRO),1)
FUN_MIRO_OBJECTS   = fun_miro/miro.io.o
MIROFLAGS          = -I$(MIRO_ROOT)/src -I$(MIRO_ROOT)/idl -I$(ACE_ROOT) -I$(TAO_ROOT) -I$(TAO_ROOT)/orbsvcs -I$(TAO_ROOT)/orbsvcs/orbsvcs \
                     -I$(QTDIR)/include -DQT_GENUINE_STR
                     # -I$(NNSIM_ROOT)/src
MIROLIBS           = -L$(MIRO_ROOT)/lib -lmiro -lmiroSvc -lmiroIDL
MIROEXTRA          = libNNSIM.a
endif

ifeq ($(USE_MIRROR),1)
FUN_MIRROR_OBJECTS = fun_miro/mirror.io.o
MIRRORFLAGS        =
MIRRORLIBS         = ../xtra/FruitFromUlrich/.obj/ObjectVectorC.o ../xtra/FruitFromUlrich/.obj/ObjectVectorI.o ../xtra/FruitFromUlrich/.obj/ObjectVectorS.o
endif

CFLAGS = $(MIROFLAGS) -DGAZEBO=$(USE_GAZEBO) -DMIRO=$(USE_MIRO)
# -c not needed


$(TARGET): $(PARSER_OBJECTS) $(KERNEL_OBJECTS) $(FUN_BASIC_OBJECTS) $(FUN_SPECIAL_OBJECTS) $(FUN_DATA_OBJECTS) $(FUN_GAZEBO_OBJECTS) $(FUN_MIRO_OBJECTS) $(FUN_MIRROR_OBJECTS)      # $(BTTVOBJECTS) $(MIROEXTRA)
	$(CC) $(LDLIBS) -DBTTV=$(USE_BTTV) -DMIRO=$(USE_MIRO) -o $(TARGET) \
	$(PARSER_OBJECTS) $(KERNEL_OBJECTS) $(FUN_BASIC_OBJECTS) $(FUN_SPECIAL_OBJECTS) $(FUN_DATA_OBJECTS) $(FUN_GAZEBO_OBJECTS) $(FUN_MIRO_OBJECTS) $(FUN_MIRROR_OBJECTS) $(GAZEBOLIBS) $(MIROLIBS) $(MIRRORLIBS) -lfl
# $(BTTVOBJECTS) $(MIROOBJECTS)

# if lex instead of flex then -ll -ly instead of -lfl


# also makes "prae.exe"
parser/y.tab.o:   parser/r.yacc.c parser/r.lex kernel/coco.h kernel/series.h
	(cd parser; make; cd ../)

libNNSIM.a: NNsimObject.o
	ar -r libNNSIM.a  $(BASICOBJECTS) $(MIROOBJECTS) $(BTTVOBJECTS) NNsimObject.o

NNsimObject.o: NNsimObject.cpp NNsimObject.h coco.h series.h data.h parser/r.yacc.h parser/y.tab.h relax.h vehicle.h observe.h utils.h
	g++ $(CFLAGS) -o NNsimObject.o NNsimObject.cpp

clean:
	\rm -f */*.o ../tools/fft/*.o */*~ $(TARGET) $(MIROEXTRA) parser/y.output parser/y.tab.h parser/y.tab.c parser/y.code.c parser/lex.yy.c ../prae.exe




# kernel/coco.o: kernel/coco.c  parser/r.yacc.h kernel/coco.h kernel/series.h fun_data/data.h kernel/relax.h kernel/vehicle.h kernel/parameters.h fun_basic/weight.h fun_basic/observe.h
# 	$(CC) $(CFLAGS) -o kernel/coco.o kernel/coco.c
# kernel/utils.o: kernel/utils.c
# 	$(CC) $(CFLAGS) -o kernel/utils.o kernel/utils.c
# data.images.o: data.images.c   coco.h series.h data.h
# 	$(CC) $(CFLAGS) -o data.images.o data.images.c
# data.digits.o: data.digits.c  coco.h series.h data.h
# 	$(CC) $(CFLAGS) -o data.digits.o data.digits.c
# data.made.o: data.made.c  coco.h series.h data.h
# 	$(CC) $(CFLAGS) -o data.made.o data.made.c
# data.bttv.o: data.bttv.c   coco.h series.h data.h
# 	$(CC) $(CFLAGS) -DBTTV=$(USE_BTTV) -o data.bttv.o data.bttv.c
# miro.io.o: miro.io.c  miro.io.h coco.h
# 	$(CPP) $(CFLAGS) $(MIROFLAGS) -o miro.io.o miro.io.c
# mirror.io.o: mirror.io.c  mirror.io.h coco.h
# 	$(CPP) $(CFLAGS) $(MIROFLAGS) $(MIRRORFLAGS) -o mirror.io.o mirror.io.c
# data.motor.o: data.motor.c  coco.h series.h data.h
# 	$(CC) $(CFLAGS) -o data.motor.o data.motor.c
# kernel/vehicle.o: kernel/vehicle.c kernel/coco.h kernel/series.h kernel/utils.h fun_basic/single.h fun_basic/local.h fun_basic/weight.h fun_basic/total.h fun_basic/observe.h \
#                   fun_special/reinforcement.h fun_data/data.h fun_miro/miro.io.h
# 	$(CC) $(CFLAGS) -DBTTV=$(USE_BTTV) -DMIRO=$(USE_MIRO) -o kernel/vehicle.o kernel/vehicle.c
# kernel/relax.o: kernel/relax.c  kernel/coco.h kernel/series.h kernel/vehicle.h
# 	$(CC) $(CFLAGS) -o kernel/relax.o kernel/relax.c
# single.o: single.c  coco.h series.h local.h
# 	$(CC) $(CFLAGS) -o single.o single.c
# local.o: local.c  coco.h
# 	$(CC) $(CFLAGS) -o local.o local.c
# fun_weight_list/weight_list.o: fun_weight_list/weight_list.c fun_weight_list/weight_list.h kernel/coco.h kernel/series.h
# 	$(CC) $(CFLAGS) -o fun_weight_list/weight_list.o fun_weight_list/weight_list.c
# total.o: total.c coco.h series.h utils.h
# 	$(CC) $(CFLAGS) -o total.o total.c
# observe.o: observe.c coco.h series.h
# 	$(CC) $(CFLAGS) -o observe.o observe.c
