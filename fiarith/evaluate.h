/*****************************************************************************/
/* evaluate.h — fiarith 表达式求值器头文件                                 */
/*****************************************************************************/

#ifndef __FIARITH_EVALUATE_H_INCLUDED
#define __FIARITH_EVALUATE_H_INCLUDED   1

#include "image.h"
#include "mask.h"
#include <psn/psn.h>

/* PSN 符号表 */
extern psnsym psn_img_op[];
extern psnsym psn_img_fn[];
extern psnprop psn_img_prop[];
extern psnsym psn_subimg_var[];

/* 操作数 */
typedef struct {
    image   *img;
} operand;

/* 求值上下文 */
typedef struct {
    int     sx, sy;
    operand *ops;
    int     nop;
    psn     **udfs;
    int     nudf;
    char    **mask;
} evaldata;

/* 栈元素 */
typedef struct {
    image   *img;
    double  val;
} imgstack;

/* 图像内存分配 */
image * image_new_empty(int sx, int sy);
image * image_new_constant(int sx, int sy, double d);

/* 核心求值器 */
image * evaluate(evaldata *ed, psn *pseq);

#endif
