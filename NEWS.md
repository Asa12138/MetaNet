# MetaNet v0.2.5 Notes

## Others

- Reconstruct the 'coors' object, now it is a dataframe rather than list. <2025-04-07, Mon>

# MetaNet v0.2.2 Notes

## Added

- added `g_layout_multi_layer` function <2025-04-04, Fri>
- added `transform_coors` function <2025-04-03, Thu>
- added `cyjs` format for function `c_net_load` <2024-11-19, Tue>

# MetaNet v0.2.1 Notes

## Added

- added `clean_multi_edge_metanet` and modified the `skeleton_plot`, now `c_net_plot` can plot skeleton network perfectly. <2024-06-05, Wed>
- added `remove_negative` in `extract_sample_net` <2024-05-23, Thu>
- added `cal_KLD` <2024-05-23, Thu>
- added plot pie function in `c_net_plot` <2024-04-12, Fri>
- added `metanet_shapes()`, including diamond, triangle1-2 <2024-04-12, Fri>
- v0.2.1, big change <2024-04-11, Thu>

# MetaNet v0.1.3 Notes

## Added

- added `c_net_union` <2024-03-29, Fri>
- added `params_list` argument for `c_net_plot`, easy to repeat calls to set default parameters <2024-03-28, Thu>

## Fixed

- fixed `zp_analyse()` as the Zi or Pi maybe NA <2024-04-09, Tue>

# MetaNet v0.1.2 Notes

## Fixed

- fixed all issues on CRAN, submitted again <2024-03-21, Thu>

# MetaNet v0.1.1 Notes

## Added

- added treemap, backbone, and stress mode for `g_layout_nice()` <2024-02-02, Fri>
- added layouts from `ggraph` to `c_net_lay()` <2024-02-02, Fri>
- added `as_polycircle()` and `as_circle_tree()` for the graph layout. <2024-02-02, Fri>
- added `module_label()` for `c_net_plot()`. <2024-01-19, Fri>
- refactored `c_net_plot()`, it's too long and complex <2024-01-31, Fri>

# MetaNet v0.1.0 Notes

## Fixed

- fixed R CMD check and Bioc check <2024-01-10>
- submitted to CRAN <2024-01-10>
- available on CRAN <https://CRAN.R-project.org/package=MetaNet> <2024-01-12>

