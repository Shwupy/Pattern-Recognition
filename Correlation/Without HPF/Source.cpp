#include <stdio.h>
#include <math.h>
#include "io_bmp.h"
#include "image_comps.h"

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
    int r, c;

    // First extend upwards, copying values of top row of image for each border row above
    float* first_line = buf;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            first_line[-r * stride + c] = first_line[c];

    // Now extend downwards, copying values of bottom row of image for each border row below
    float* last_line = buf + (height - 1) * stride;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            last_line[r * stride + c] = last_line[c];

    // Now extend all rows to the left and to the right by the value in the left/right_edge[0]
    float* left_edge = buf - border * stride;
    float* right_edge = left_edge + width - 1;
    for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
        for (c = 1; c <= border; c++)
        {
            left_edge[-c] = left_edge[0];
            right_edge[c] = right_edge[0];
        }
}

/*****************************************************************************/
/*                              External Functions                           */
/*****************************************************************************/
float J_COR(float* search, int search_width, float* pattern, int pat_height, int pat_width) {
    //"search" is the location within the search image we are checking
    //"pattern" is the top left/beginning of pattern image we're using
    float sum = 0;
    float* tempp = pattern;
    float* temps = search;
    for (int r = 0; r < pat_height; r++) {
        for (int c = 0; c < pat_width; c++) {
            tempp = pattern + r * pat_width + c;    //pat[m][n]
            temps = search + r * search_width + c;  //search[x+m][y+n]
            float multiply = (*tempp) * (*temps);
            sum = sum + multiply;
            //printf("At coord: [%d] [%d], sqr_diff = %f\n", r, c, sqr_diff);
            //printf("Sum = %f\n", sum);
        }
    }
    return sum;
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char* argv[])
{
    if (argc != 3) //exe name, 2 input files
    {
        fprintf(stderr, "Usage: %s <search bmp file> <pattern bmp file>\n", argv[0]);
        return -1;
    }

    int err_code = 0;
    try {
        // Read the 1st input image
        bmp_in in;
        if ((err_code = bmp_in__open(&in, argv[1])) != 0)
            throw err_code;

        int width = in.cols, height = in.rows;
        int n, num_comps = in.num_components;
        my_image_comp* input_comps = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps[n].init(height, width, 0); // Leave a border of 4

        int r; // Declare row index
        io_byte* line = new io_byte[width * num_comps];
        for (r = height - 1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
            if ((err_code = bmp_in__get_line(&in, line)) != 0)
                throw err_code;
            for (n = 0; n < num_comps; n++)
            {
                io_byte* src = line + n; // Points to first sample of component n
                float* dst = input_comps[n].buf + r * input_comps[n].stride;
                for (int c = 0; c < width; c++, src += num_comps)
                    dst[c] = (float)*src; // The cast to type "float" is not
                          // strictly required here, since bytes can always be
                          // converted to floats without any loss of information.
            }
        }
        bmp_in__close(&in);

        // Read the 2nd input image
        bmp_in in1;
        if ((err_code = bmp_in__open(&in1, argv[2])) != 0)
            throw err_code;

        int width1 = in1.cols, height1 = in1.rows;
        int n1, num_comps1 = in1.num_components;
        my_image_comp* input_comps1 = new my_image_comp[num_comps1];
        for (n1 = 0; n1 < num_comps1; n1++)
            input_comps1[n1].init(height1, width1, 0); // Leave a border of 4

        int r1; // Declare row index
        io_byte* line1 = new io_byte[width1 * num_comps1];
        for (r1 = height1 - 1; r1 >= 0; r1--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
            if ((err_code = bmp_in__get_line(&in1, line1)) != 0)
                throw err_code;
            for (n1 = 0; n1 < num_comps1; n1++)
            {
                io_byte* src1 = line1 + n1; // Points to first sample of component n
                float* dst1 = input_comps1[n1].buf + r1 * input_comps1[n1].stride;
                for (int c = 0; c < width1; c++, src1 += num_comps1)
                    dst1[c] = (float)*src1; // The cast to type "float" is not
                          // strictly required here, since bytes can always be
                          // converted to floats without any loss of information.
            }
        }
        bmp_in__close(&in1);

        if (height < height1 || width < width1) {
            printf("Search image must be larger than Pattern image\n");
            return -1;
        }
        printf("Search image has size: height %d by width %d\n", height, width);
        printf("Pattern image has size: height %d by width %d\n", height1, width1);
        float min = 0;
        int pix_x = 0;
        int pix_y = 0;
        for (int r = 0; r < height - height1; r++) {
            for (int c = 0; c < width - width1; c++) {
                float* search_coord = input_comps[0].buf + r * width + c;
                float result = J_COR(search_coord, width, input_comps1[0].buf, height1, width1);
                //printf("result = %f\n", result);
                if (result > min) {
                    min = result;
                    //printf("result = %f\n", result);
                    pix_y = r;
                    pix_x = c;
                }
            }
        }
        printf("matching pixel at (x,y) =  (%d,%d)\n", pix_x, pix_y);
        delete[] line;
        delete[] input_comps;
        delete[] input_comps1;
    }
    catch (int exc) {
        if (exc == IO_ERR_NO_FILE)
            fprintf(stderr, "Cannot open supplied input or output file.\n");
        else if (exc == IO_ERR_FILE_HEADER)
            fprintf(stderr, "Error encountered while parsing BMP file header.\n");
        else if (exc == IO_ERR_UNSUPPORTED)
            fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
        else if (exc == IO_ERR_FILE_TRUNC)
            fprintf(stderr, "Input or output file truncated unexpectedly.\n");
        else if (exc == IO_ERR_FILE_NOT_OPEN)
            fprintf(stderr, "Trying to access a file which is not open!(?)\n");
        return -1;
    }
    return 0;
}
