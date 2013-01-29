// Contents of DLLDefines.h
#ifndef _CAPGMath_DLLDEFINES_H_
#define _CAPGMath_DLLDEFINES_H_

/* Cmake will define CAPGMath_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define CAPGMath_EXPORTS when
building a DLL on windows.
*/
// We are using the Visual Studio Compiler and building Shared libraries

#if defined (_WIN32) 
  #if defined(CAPGMath_EXPORTS)
    #define  _CAPGExport __declspec(dllexport)
  #else
    #define  _CAPGExport __declspec(dllimport)
  #endif /* CAPGMath_EXPORTS */
#else /* defined (_WIN32) */
 #define _CAPGExport
#endif

#endif /* _CAPGMath_DLLDEFINES_H_ */

