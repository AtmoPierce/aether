#include <math.h>
#include <stdint.h>

double fpx_cos_x87(double x) {
    double y;
    __asm__ volatile (
        "fldl  %[x]\n\t"
        "fcos\n\t"
        "fstpl %[y]\n\t"
        : [y] "=m"(y)
        : [x] "m"(x)
        : "st"
    );
    return y;
}
