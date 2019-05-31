LEX  = flex -I
YACC = yacc -d -v        # byacc on PC;  option -t on computers other than gort

                         # note about file r.lex:
                         # for gcc use input()
                         # for g++ use yyinput()

all:    lex yacc prae.exe
	g++ -c lex.yy.c y.tab.c

lex:    r.lex
	${LEX} r.lex

yacc:   r.yacc.c
	${YACC} r.yacc.c

prae.exe: prae.c
	g++ -Wall prae.c -o ../../prae.exe

clean:
	\rm -f *.o *~ y.output y.tab.h y.tab.c lex.yy.c
