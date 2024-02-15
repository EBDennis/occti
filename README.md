# occti R package
R package to produce species' occupancy trends and indices using the unmarked package. Occupancy models are applied separately to each year which are run in reverse order, with starting values informed by the preceding year.

The approach is based upon: Dennis, E.B., Morgan, B.J.T., Freeman, S.N., Ridout, M.S., Brereton, T.M., Fox, R., Powney, G.D. and Roy, D.B. (2017) Efficient occupancy model-fitting for extensive citizen-science data. _Plos One_, 12(3): e0174433. https://doi.org/10.1371/journal.pone.0174433

Calculation of standard errors using the Delta method follows Dennis, E.B., Brereton, T.M., Morgan, B.J.T., Fox, R., Shortall, C.R., Prescott, T. and Foster, S. (2019) Trends and indicators for quantifying moth abundance and occupancy in Scotland. _Journal of Insect Conservation_, 23, 369–380. https://doi.org/10.1007/s10841-019-00135-z
