%{
   #include <stdlib.h>
   #include <stdio.h>
   #include <string.h>
   #include "y.tab.h"

   extern int lineno;
   extern void warning(char *s, char *t);
%}

ws              [ \t]+
opencomment     #.*$
id              [a-zA-Z\_\.-/][a-zA-Z\_\.0-9/]*
number          (-?[0-9]+)|-?([0-9]*\.[0-9]+)
gi_             \{
_ig             \}
ri_             \(
_ir             \)
quote           \"
plus            \+
mult            \*
colon           \:
newline         \n
sep             ;
komma           ,

%%

(global|GLOBAL|Global)    { return GLOBAL; }
(all|ALL|All)             { return ALL; }
(area|AREA|Area)          { return AREAS; }
(series|SERIES|Series)    { return SER; }
(sw|SW|sweep|SWEEP|Sweep) { return SWE; }
(if|IF|If)                { return IFWORD; }
"="                       { return EQ; }
"%"                       { return MOD; }
"<"                       { return SM; }
">"                       { return LG; }
(iter)                    { return ITERWORD; }
"/*"                      { int c1 = 0, c2 = yyinput(); /*use input() with gcc; use yyinput() with g++*/
                            for (;;) {
                                if  (c2 == EOF)
                                    break;
                                if  ((c1 == '*') && (c2 == '/'))
                                    break;
                                c1 = c2;
                                c2 = yyinput();
                            }
                          }
"["                       { int c1 = yyinput();
                            for (;;) {
                                if  (c1 == EOF)
                                    break;
                                if  (c1 == ']')
                                    break;
                                c1 = yyinput();
                            }
                          }
{ws}                      ;
{opencomment}             { printf("COMMENT: %s\n",yytext); }
{sep}                     { return SEP; }
{komma}                   { return KOMMA; }
{gi_}                     { return GI_; }
{_ig}                     { return _IG; }
{ri_}                     { return RI_; }
{_ir}                     { return _IR; }
{quote}                   { return QUOTE; }
{plus}                    { return PLUS; }
{mult}                    { return MULT; }
{colon}                   { return COLON; }
{newline}                 { lineno++; }
{number}                  { yylval.string = strdup(yytext);
			    return NUMBER; 
			  }
{id}                      { yylval.string = strdup(yytext);
			    return ID; 
			  }
.                         { printf("\nnot in lex table:%c\n", yytext[0]);
                            return yytext[0];
                          }

%%
