@ECHO OFF

REM set ProblemName     = %1 
REM set ProblemPath     = %2
REM set ProblemTypePath = %3

REM OutputFile: %2\%1.info

del %2\%1.info
del %2\%1.flavia.res
del %2\%1.post.res
del %2\%1.post.msh

%3\execs\win\ERMES_20.0.4.exe %1 > %2\%1.info