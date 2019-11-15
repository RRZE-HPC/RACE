find . ! \( -path ./matrices/table.tex -prune -o -path \
./results/alpha_table_generator/table.tex -prune -o -path \
./level_construction/FDM_2d_7pt_non_perm.png -prune -o -path \
./permutation/FDM_2d_7pt_perm.png -prune -o -name "*.pdf" \) -delete
