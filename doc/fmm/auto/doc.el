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
    "doi")
   (LaTeX-add-labels
    "eq:HJ"
    "eq:eikonal"
    "eq:eikonal_discr"
    "eq:osher"
    "eq:eikonal_discr_2")
   (LaTeX-add-bibitems
    "rickett99"
    "rouy92"
    "sethian96"
    "sethian99"
    "sethian03"
    "thomsen86"
    "tsai03"
    "tsitsiklis95")))

