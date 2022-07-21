#include "sparsemat.h"
#include "multicoloring.h"
#include <vector>
#include <map>
#include "config_eg.h"
//uses METIS for partitioning
#ifdef RACE_HAVE_METIS
    #include "metis.h"
#endif
#ifdef RACE_HAVE_COLPACK
    #include "ColPackHeaders.h"
#else
//    #warning "Coloring needs COLPACK. Please compile and link with COLPACK if you use coloring"
#endif
#include <limits>
#include "timer.h"

multicoloring::multicoloring(sparsemat *mat_, int colorDist_, int* initPerm_, int* initInvPerm_):mat(mat_), colorDist(colorDist_), initPerm(initPerm_), initInvPerm(initInvPerm_), colorBlockSize_out(-1), ncolors_out(-1), colorPtr_out(NULL), partPtr_out(NULL), perm_out(NULL), invPerm_out(NULL)
{
    if(colorDist != 1 && colorDist != 2)
    {
        printf("Dist %d not supported, falling back to dist-2\n",colorDist);
        colorDist = 2;
    }
}

bool multicoloring::doMC()
{
#ifdef RACE_HAVE_COLPACK
    ColPack::GraphColoringInterface *GC=new ColPack::GraphColoringInterface(-1);
    int nnz = mat->nnz;
    int nrows = mat->nrows;
    int* rptlocal = (int*) malloc(sizeof(int)*(nrows+1));
    int* collocal = (int*) malloc(sizeof(int)*nnz);
    rptlocal[0] = 0;
    int perm_nnz = 0;
    for(int row=0; row<nrows; ++row)
    {
        int permRow = row;
        if(initPerm)
        {
            permRow = initPerm[row];
        }
        for(int idx=mat->rowPtr[permRow]; idx<mat->rowPtr[permRow+1]; ++idx)
        {
            int permCol = mat->col[idx];
            if(initInvPerm)
            {
                permCol = initInvPerm[mat->col[idx]];
            }
            collocal[perm_nnz] = permCol;
            ++perm_nnz;
        }
        rptlocal[row+1] = perm_nnz;
    }

    uint32_t** adolc = new uint32_t*[nrows];
    uint32_t* adolc_data = new uint32_t[nnz+nrows];
    int64_t pos=0;
    for (int i=0;i<nrows;i++)
    {
        adolc[i]=&(adolc_data[pos]);
        adolc_data[pos++]=rptlocal[i+1]-rptlocal[i];
        for (int j=rptlocal[i];j<rptlocal[i+1];j++)
        {
            adolc_data[pos++]=collocal[j];
        }
    }
    GC->BuildGraphFromRowCompressedFormat(adolc, nrows);
    if(colorDist == 1)
    {
        int ret = GC->Coloring();
        if(ret==false)
        {
            printf("Error in COLPACK coloring\n");
        }
    }
    else
    {
        int ret = GC->Coloring("NATURAL", "DISTANCE_TWO");
        if(ret==false)
        {
            printf("Error in COLPACK coloring\n");
        }

        //COLPACK_CALL_GOTO(GC->DistanceTwoColoring(),err,ret);
        /*
           if (GC->CheckDistanceTwoColoring(2)) {
           GHOST_ERROR_LOG("Error in coloring!");
           ret = GHOST_ERR_COLPACK;
           goto err;
           }*/
    }
    int ncolors = GC->GetVertexColorCount();
    printf("No. of colors = %d\n", ncolors);
    //colorPtr rowPtr of parts for each color
    int* colorPtr = (int*)malloc(sizeof(int)*(ncolors+1));
    //stores the color of each part
    std::vector<int>* colvec = NULL;
    colvec = GC->GetVertexColorsPtr();

    for (int i=0;i<ncolors+1;i++) {
        colorPtr[i] = 0;
    }

    for (int i=0;i<nrows;++i) {
        colorPtr[(*colvec)[i]+1]++;
    }

    for (int i=1;i<ncolors+1;i++) {
        colorPtr[i] += colorPtr[i-1];
    }

    int* curcol = (int*)malloc(sizeof(int)*(ncolors));
    memset(curcol,0,ncolors*sizeof(int));

    /*permute according to color*/
    int *finalPerm = (int*)malloc(sizeof(int)*nrows);
    int *finalInvPerm = (int*)malloc(sizeof(int)*nrows);
    for(int i=0; i<nrows; ++i)
    {
        finalPerm[i]=0;
        finalInvPerm[i]=0;
    }
    for(int row=0; row<nrows; ++row)
    {
        //find rows in parts
        int permCurrRow = row;
        if(initPerm)
        {
            permCurrRow = initPerm[row];
        }
        finalInvPerm[permCurrRow] = curcol[(*colvec)[row]] + colorPtr[(*colvec)[row]];
        curcol[(*colvec)[row]]++;
    }
    for (int i=0; i<nrows; i++) {
        finalPerm[finalInvPerm[i]] = i;
    }

    //transfer data structure to output
    colorBlockSize_out = 1;
    ncolors_out = ncolors;
    if(colorPtr_out)
    {
        free(colorPtr_out);
    }
    colorPtr_out = colorPtr;
    if(perm_out)
    {
        free(perm_out);
    }
    perm_out = finalPerm;
    if(invPerm_out)
    {
        free(invPerm_out);
    }
    invPerm_out = finalInvPerm;

    delete [] adolc_data;
    delete [] adolc;
    delete GC;
    free(rptlocal);
    free(collocal);
    free(curcol);
    return 0;
#else
    printf("ColPack not available to perform MC. Please compile and link with ColPack\n");
    return -1;
#endif

}

bool multicoloring::doABMC()
{
#if (defined RACE_HAVE_METIS && defined RACE_HAVE_COLPACK)
    char *blockSize_env = (getenv("ABMC_BLOCKSIZE"));

    int blocksize = 64;
    if(blockSize_env)
    {
        blocksize= atoi(blockSize_env);
    }

    //std::vector<int> block_sizes_default{8,32,64,128};
    std::vector<int> block_sizes_default{64,128};
    std::vector<int> block_sizes;
    int opt_block;

    if(blockSize_env)
    {
        block_sizes.push_back(blocksize);
        opt_block = blocksize;
    }
    else
    {
        block_sizes = block_sizes_default;
    }

    double min_time = std::numeric_limits<double>::max();
    for(int blockIdx=0; blockIdx<(int)block_sizes.size(); ++blockIdx)
    {
        blocksize = block_sizes[blockIdx];
        printf("blocksize = %d\n", blocksize);
        ColPack::GraphColoringInterface *GC=new ColPack::GraphColoringInterface(-1);
        int nnz = mat->nnz;
        int nrows = mat->nrows;
        int* rptlocal = (int*) malloc(sizeof(int)*(nrows+1));
        int* collocal = (int*) malloc(sizeof(int)*nnz);
        rptlocal[0] = 0;
        int perm_nnz = 0;
        for(int row=0; row<nrows; ++row)
        {
            int permRow = row;
            if(initPerm)
            {
                permRow = initPerm[row];
            }
            for(int idx=mat->rowPtr[permRow]; idx<mat->rowPtr[permRow+1]; ++idx)
            {
                int permCol = mat->col[idx];
                if(initInvPerm)
                {
                    permCol = initInvPerm[mat->col[idx]];
                }
                collocal[perm_nnz] = permCol;
                ++perm_nnz;
            }
            rptlocal[row+1] = perm_nnz;
        }

        //partition using METIS
        int ncon = 1;
        int nparts = (int)(nrows/(double)blocksize);
        int objval;
        int *part = (int*) malloc(sizeof(int)*nrows);

        printf("partitioning graph to %d parts\n", nparts);
        int metis_ret = METIS_PartGraphKway(&nrows, &ncon, rptlocal, collocal, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);
        if(metis_ret == METIS_OK)
        {
            printf("successfully partitioned graph to nparts=%d\n", nparts);
        }
        else
        {
            printf("Error in ABMC partitioning\n");
        }

        std::vector<std::vector<int>> partRow(nparts);
        int *partPtr = (int*)malloc(sizeof(int)*(nparts+1));
        for (int i=0;i<nparts+1;i++) {
            partPtr[i] = 0;
        }
        for (int i=0;i<nrows;i++) {
            partRow[part[i]].push_back(i);
            partPtr[part[i]+1]++;
        }
        for (int i=1;i<nparts+1;i++) {
            partPtr[i] += partPtr[i-1];
        }
        //#now make a new matrix from the partitions
        int nrows_b = nparts;
        int *rowPtr_b = (int*)malloc(sizeof(int)*(nrows_b+1));
        for(int i=0; i<(nrows_b+1); ++i)
        {
            rowPtr_b[i] = 0;
        }

        //unique_col stores partners parts to which the current part is connected
        std::vector<std::map<int,int>> unique_col(nrows_b);
        //for each part I need to know to which part it is connected
        for(int row=0; row<nrows; ++row)
        {
            int currPart_id = part[row];
            for(int idx=rptlocal[row]; idx<rptlocal[row+1]; ++idx)
            {
                int currCol_idx = collocal[idx];
                //find to which part does this col belong
                int partnerPart = part[currCol_idx];
                unique_col[currPart_id][partnerPart]=1;//I just need to create unique entries; therefore map
            }
        }

        std::vector<int> col_b;
        //count nnzr
        for(int part_idx=0; part_idx<nrows_b; ++part_idx)
        {
            int nnzr_b = unique_col[part_idx].size();
            rowPtr_b[part_idx+1] = rowPtr_b[part_idx] + nnzr_b;
            for(auto& map_el : unique_col[part_idx])
            {
                col_b.push_back(map_el.first);
            }
        }

        //now color this graph
        int nnz_b = col_b.size();
        uint32_t** adolc = new uint32_t*[nrows_b];
        uint32_t* adolc_data = new uint32_t[nnz_b+nrows_b];
        int64_t pos=0;
        for (int i=0;i<nrows_b;i++)
        {
            adolc[i]=&(adolc_data[pos]);
            adolc_data[pos++]=(rowPtr_b[i+1]-rowPtr_b[i]);
            for (int j=rowPtr_b[i];j<rowPtr_b[i+1];j++)
            {
                adolc_data[pos++]=col_b[j];
            }
        }

        GC->BuildGraphFromRowCompressedFormat(adolc, nrows_b);
        if(colorDist == 1)
        {
            int ret = GC->Coloring();
            if(ret==false)
            {
                printf("Error in COLPACK coloring\n");
            }
        }
        else
        {
            int ret = GC->Coloring("NATURAL", "DISTANCE_TWO");
            if(ret==false)
            {
                printf("Error in COLPACK coloring\n");
            }

            //COLPACK_CALL_GOTO(GC->DistanceTwoColoring(),err,ret);
            /*
               if (GC->CheckDistanceTwoColoring(2)) {
               GHOST_ERROR_LOG("Error in coloring!");
               ret = GHOST_ERR_COLPACK;
               goto err;
               }*/
        }

        int ncolors = GC->GetVertexColorCount();
        printf("No. of colors = %d\n", ncolors);

        //stores rowPtr of actual row for each color
        int* color_part_ptr = (int*)malloc(sizeof(int)*(ncolors+1));
        //colorPtr rowPtr of parts for each color
        int* colorPtr = (int*)malloc(sizeof(int)*(ncolors+1));
        //stores the color of each part
        std::vector<int>* colvec = NULL;
        colvec = GC->GetVertexColorsPtr();

        for (int i=0;i<ncolors+1;i++) {
            colorPtr[i] = 0;
            color_part_ptr[i] = 0;
        }

        for (int i=0;i<nrows_b;++i) {
            colorPtr[(*colvec)[i]+1]++;
            //the size of current part
            int partSize = partPtr[i+1] - partPtr[i];
            //How much part contain this color
            color_part_ptr[(*colvec)[i]+1]+=partSize;
        }

        for (int i=1;i<ncolors+1;i++) {
            colorPtr[i] += colorPtr[i-1];
            color_part_ptr[i] += color_part_ptr[i-1];
        }

        int* curcol = (int*)malloc(sizeof(int)*(ncolors));
        memset(curcol,0,ncolors*sizeof(int));
        std::vector<std::vector<int>> color_part_map(ncolors);

        int *finalPerm = (int*)malloc(sizeof(int)*nrows);
        int *finalInvPerm = (int*)malloc(sizeof(int)*nrows);
        for(int i=0; i<nrows; ++i)
        {
            finalPerm[i]=0;
            finalInvPerm[i]=0;
        }
        for (int partIdx=0; partIdx<nrows_b; partIdx++)
        {
            color_part_map[(*colvec)[partIdx]].push_back(partIdx);
            int partSize = partPtr[partIdx+1]-partPtr[partIdx];
            for(int rowIdx=0; rowIdx<partSize; ++rowIdx)
            {
                //find rows in parts
                int currRow = partRow[partIdx][rowIdx];
                int permCurrRow = currRow;
                if(initPerm)
                {
                    permCurrRow = initPerm[currRow];
                }
                finalInvPerm[permCurrRow] = curcol[(*colvec)[partIdx]] + color_part_ptr[(*colvec)[partIdx]];
                curcol[(*colvec)[partIdx]]++;
            }
        }
        for (int i=0; i<nrows; i++) {
            finalPerm[finalInvPerm[i]] = i;
        }

        //permute the partPtr according to color
        int* permPartPtr = (int*)malloc(sizeof(int)*(nrows_b+1));
        permPartPtr[0] =0;
        int ctr=0;
        //now permute the partPtr
        for(int chrom=0; chrom<ncolors; ++chrom)
        {
            for(int k=0; k<(int)color_part_map[chrom].size(); ++k)
            {
                int currPartIdx = color_part_map[chrom][k];
                permPartPtr[ctr+1] = partPtr[currPartIdx+1]-partPtr[currPartIdx];
                ++ctr;
            }
        }

        for(int partIdx=0; partIdx<nrows_b; ++partIdx)
        {
            permPartPtr[partIdx+1] += permPartPtr[partIdx];
        }


        //create matrix  and test performance
        int *new_rptlocal = (int*) malloc(sizeof(int)*(nrows+1));
        int *new_collocal = (int*) malloc(sizeof(int)*nnz);
        double *new_vallocal = (double*) malloc(sizeof(double)*nnz);
        double *b = (double*) malloc(sizeof(double)*nrows);
        double *x = (double*) malloc(sizeof(double)*nrows);

        new_rptlocal[0] = 0;
        int nnzlocal = 0;

        for (int i=0; i<nrows; i++) {
            b[i] = 0;
            x[i] = 0.001;
            new_rptlocal[i+1] = new_rptlocal[i];

            int orig_row = i;
            if (finalPerm) {
                orig_row = finalPerm[i];
                // orig_row = ctx->row_map->loc_perm[i];
            }
            int * col = &mat->col[mat->rowPtr[orig_row]];
            int orig_row_len = mat->rowPtr[orig_row+1]-mat->rowPtr[orig_row];

            for(int j=0; j<orig_row_len; ++j) {
                if(finalPerm)
                {
                    new_collocal[nnzlocal] = finalInvPerm[col[j]];
                    new_vallocal[nnzlocal] = 0.1;
                }
                else
                {
                    new_collocal[nnzlocal] = col[j];
                    new_vallocal[nnzlocal] = 0.1;
                }
                nnzlocal++;
                new_rptlocal[i+1]++;
            }
        }

        double start_time, end_time;

        INIT_TIMER(spmv_abmc_test);
        START_TIMER(spmv_abmc_test);
        double curTime = 0;
        int nIter = 0;
        while(curTime<0.5)
        {

            //do SpMV and check
#pragma omp parallel for schedule(static)
            for (int row=0; row<nrows; ++row) {
                double temp = 0;
                int idx = new_rptlocal[row];
#pragma simd reduction(+:temp)
                for (int j=new_rptlocal[row]; j<new_rptlocal[row+1]; j++) {
                    temp = temp + new_vallocal[j] * x[new_collocal[j]];
                }
                b[row]=temp;
            }

            STOP_TIMER(spmv_abmc_test);
            curTime = GET_TIMER(spmv_abmc_test);
            //swap pointers
            double* tmp = x;
            x = b;
            b = tmp;
            ++nIter;
        }

        double tot_time = curTime/(double)nIter;
        printf("cur time = %f s\n", tot_time);
        if(tot_time<min_time)
        {
            min_time = tot_time;
            opt_block = blocksize;

            //store relevant data structures
            colorBlockSize_out = opt_block;
            ncolors_out = ncolors;
            if(colorPtr_out)
            {
                free(colorPtr_out);
            }
            colorPtr_out = colorPtr;
            if(partPtr_out)
            {
                free(partPtr_out);
            }
            partPtr_out = permPartPtr;
            if(perm_out)
            {
                free(perm_out);
            }
            perm_out=finalPerm;
            if(invPerm_out)
            {
                free(invPerm_out);
            }
            invPerm_out=finalInvPerm;
        }
        else
        {
            free(colorPtr);
            colorPtr=NULL;
            free(permPartPtr);
            permPartPtr=NULL;
            free(finalPerm);
            finalPerm=NULL;
            free(finalInvPerm);
            finalInvPerm=NULL;
        }

        free(new_rptlocal);
        free(new_collocal);
        free(new_vallocal);
        free(b);
        free(x);

        delete [] adolc_data;
        delete [] adolc;
        delete GC;
        // free(rpt);
        free(rptlocal);
        free(collocal);
        free(part);
        free(partPtr);
        free(rowPtr_b);
        free(color_part_ptr);
        free(curcol);
    }
    printf("Chosen block size = %d\n", colorBlockSize_out);
    return 0;
#else
    printf("METIS or ColPack not available to perform ABMC. Please compile and link with ColPack and METIS\n");
    return -1;
#endif
}

