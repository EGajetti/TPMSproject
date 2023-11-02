
#ifndef VTKBOOL_EXPORT_H
#define VTKBOOL_EXPORT_H

#ifdef VTKBOOL_STATIC_DEFINE
#  define VTKBOOL_EXPORT
#  define VTKBOOL_NO_EXPORT
#else
#  ifndef VTKBOOL_EXPORT
#    ifdef vtkBool_EXPORTS
        /* We are building this library */
#      define VTKBOOL_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define VTKBOOL_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef VTKBOOL_NO_EXPORT
#    define VTKBOOL_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef VTKBOOL_DEPRECATED
#  define VTKBOOL_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef VTKBOOL_DEPRECATED_EXPORT
#  define VTKBOOL_DEPRECATED_EXPORT VTKBOOL_EXPORT VTKBOOL_DEPRECATED
#endif

#ifndef VTKBOOL_DEPRECATED_NO_EXPORT
#  define VTKBOOL_DEPRECATED_NO_EXPORT VTKBOOL_NO_EXPORT VTKBOOL_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VTKBOOL_NO_DEPRECATED
#    define VTKBOOL_NO_DEPRECATED
#  endif
#endif

#endif /* VTKBOOL_EXPORT_H */
