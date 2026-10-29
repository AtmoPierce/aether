#include <math.h>
#include <stdint.h>

double fpx_tan_x87(double x) {
    double y;
    __asm__ volatile (
        "fldl %[x]\n\t"
        "fptan\n\t"
        "fstp %%st(0)\n\t"   // pop 1.0
        "fstpl %[y]\n\t"
        : [y] "=m"(y)
        : [x] "m"(x)
        : "st"
    );
    return y;
}
