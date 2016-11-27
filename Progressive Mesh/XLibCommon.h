#ifndef XLIB_COMMON_H_
#define XLIB_COMMON_H_

#include <iostream>
#include <string>
#include <sstream>
namespace XMeshLib {

    typedef double xreal;

    #define PI 3.1415926535897932384626

    static std::string GetExtensionName(const char filename[]) {
        std::string cstr(filename);
        std::string extName;
        int pos = cstr.find_last_of(".");
        if (pos > 0)
            extName = cstr.substr(pos + 1);
        return extName;
    }

    static std::string GenerateIndexedFileName(const char prefix[], int index, const char extName[]) {
        std::string pstr(prefix);
        std::string estr(extName);
        // std::string nstr = std::to_string(index);
        return pstr +  estr;  // nstr +
    }

    static void errorExit(const char msg[]) {
        std::cerr << msg << "\n";
        exit(0);
    }

}
#endif
