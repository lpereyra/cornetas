#ifndef COLORES_H
#define COLORES_H

extern char message[200];

#define YELLOW(a)   fprintf(stdout,"\033[01;33m%s\033[22;0m",(a))
#define DEFA(a)     fprintf(stdout,"\033[01;37m%s\033[22;0m",(a))
#define BLACK(a)    fprintf(stdout,"\033[22;30m%s\033[22;0m",(a))
#define RED(a)      fprintf(stdout,"\033[22;31m%s\033[22;0m",(a))
#define GREEN(a)    fprintf(stdout,"\033[22;32m%s\033[22;0m",(a)) 
#define BROWN(a)    fprintf(stdout,"\033[22;33m%s\033[22;0m",(a)) 
#define BLUE(a)     fprintf(stdout,"\033[22;34m%s\033[22;0m",(a)) 
#define MAGENTA(a)  fprintf(stdout,"\033[22;35m%s\033[22;0m",(a)) 
#define CYAN(a)     fprintf(stdout,"\033[22;36m%s\033[22;0m",(a)) 
#define GRAY(a)     fprintf(stdout,"\033[22;37m%s\033[22;0m",(a)) 
#define DGRAY(a)    fprintf(stdout,"\033[01;30m%s\033[22;0m",(a)) 
#define LRED(a)     fprintf(stdout,"\033[01;31m%s\033[22;0m",(a)) 
#define LGREEN(a)   fprintf(stdout,"\033[01;32m%s\033[22;0m",(a)) 
#define LBLUE(a)    fprintf(stdout,"\033[01;34m%s\033[22;0m",(a)) 
#define LMAGENTA(a) fprintf(stdout,"\033[01;35m%s\033[22;0m",(a)) 
#define LCYAN(a)    fprintf(stdout,"\033[01;36m%s\033[22;0m",(a)) 
#define WHITE(a)    fprintf(stdout,"\033[01;37m%s\033[22;0m",(a))

#endif
