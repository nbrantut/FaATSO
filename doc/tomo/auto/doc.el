(TeX-add-style-hook
 "doc"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "englbst"
    "article"
    "art10"
    "amsmath"
    "amssymb"
    "graphicx"
    "natbib")
   (TeX-add-symbols
    '("capitalize" 1)
    '("Capitalize" 1)
    '("doi" 1)
    '("natexlab" 1)
    '("mat" 1)
    "doi")
   (LaTeX-add-labels
    "eq:d=g(m)"
    "eq:Gn"
    "eq:Tij"
    "itm:interp")
   (LaTeX-add-bibitems
    "aki76"
    "tarantola05"
    "thomsen86"
    "thurber83")))

