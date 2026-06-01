---
project: hompack
license: mit
summary: hompack is a package for solving nonlinear systems of equations by homotopy methods.
src_dir: src
         example
exclude: src/lapack.f
output_dir: _site
page_dir: doc
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            iso_c_binding:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html#ISO_005fC_005fBINDING
source: true
proc_internals: true
preprocess: false
graph: true
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/HugoMVale/odrpack95
author: H. M. Vale,
        L. T. Watson,
        S. C. Billups,
        A. P. Morgan,
        M. Sosonkina,
        R. C. Melville,
        H. F. Walker
github: https://github.com/HugoMVale/
email: 57530119+HugoMVale@users.noreply.github.com
dbg: true
predocmark: >
docmark_alt: #
predocmark_alt: <
md_extensions: markdown.extensions.toc
---

`hompack` is a package for finding zeros or fixed points of nonlinear systems using globally convergent probability-one homotopy algorithms.
