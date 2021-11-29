#include <stdio.h>
#include <stdlib.h>
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
struct coord {
    int x = 0;
    int y = 0;
};
struct filt {
    float* centre;
    int length;
};
struct min_max {
    float min;
    float max;
};

//Correlation
float J_COR(float* search, int search_stride, float* pattern, int pat_height, int pat_width) {
    //"search" is the location within the search image we are checking
    //"pattern" is the top left/beginning of pattern image we're using
    float sum = 0;
    float* tempp = pattern;
    float* temps = search;
    for (int r = 0; r < pat_height; r++) {
        for (int c = 0; c < pat_width; c++) {
            tempp = pattern + r * pat_width + c;    //pat[m][n]
            temps = search + r * search_stride + c;  //search[x+m][y+n]
            float multiply = (*tempp) * (*temps);
            sum = sum + multiply;
            //printf("At coord: [%d] [%d], sqr_diff = %f\n", r, c, sqr_diff);
            //printf("Sum = %f\n", sum);
        }
    }
    return sum;
}
coord find_coord(my_image_comp* search, my_image_comp* pattern) {
    int sw = search->width;
    int sh = search->height;
    int pw = pattern->width;
    int ph = pattern->height;
    float max = 0;
    coord coords;

    //search area cut off by height and width of pattern
    for (int r = 0; r < sh-ph; r++) {
        for (int c = 0; c < sw-pw; c++) {
            float* temps = search->buf + r * sw + c; //search at position [r][c]
            float corr_val = J_COR(temps, sw, pattern->buf, ph, pw);
            //printf("%f\n", corr_val);
            if (corr_val > max) {   //Finds maximum

                max = corr_val;
                coords.x = c;
                coords.y = r;
            }
        }
    }
    return coords;
}

//Filtering
filt make_filter(int type) {
    filt filter;
    if (type == 1) {
        float* hpf1 = new float[9];
        filter.centre = hpf1 + 4;
        filter.length = 3;
        for (int i = 0; i < 9; i++) {
            if (i % 2 == 0 && i != 4) {
                hpf1[i] = -0.1406;
            }
            else if (i == 4) {
                hpf1[i] = 1.125;
            }
            else {
                hpf1[i] = -0.1406;
            }
        }
    }
    return filter;
}
float convolve(float* pattern, int patt_stride, float* filter, int H) {
    //"pattern" is the start pixel of area we want to convolve
    int length = 2 * H + 1;
    float sum = 0;

    for (int r = -H; r <= H; r++) {
        for (int c = -H; c <= H; c++) {
            sum += pattern[r * patt_stride + c] * filter[r * length + c];
        }
    }
    return sum;
}
min_max apply_hpf(float* pattern, float* out, int height, int width, int stride, float* hpf, int H) {
    //"out" needed as previous convole would take the previous convolved pixel into the next convolution
    min_max range;
    float max = 0;
    float min = 0;
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
            float* temp_coord = pattern + r * stride + c;
            *(out + r*width + c) = convolve(temp_coord, stride, hpf, H);
            //printf(" at [x][y] = [%d][%d], conv = %f\n", c, r, *(out + r * width + c));
            if (*(out + r * width + c) > max) {        //keep track of maximum value
                max = *(out + r * width + c);
            }
            if (*(out + r * width + c) < min) {        //keep track of minimum value
                min = *(out + r * width + c);
            }
        }
    }
    range.min = min;
    range.max = max;
    return range;
    //printf("min = %f, max = %f\n", min, max);
    //printf("\n After Normalising \n");
    /*
    //keep result in range between 0 and 255: add |min| then multiply 255/new_max
    float abs_min = fabs(min);
    float range = max + abs_min;   //|min| + max
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
            float* temp_coord = out + r * width + c;
            *temp_coord = make_in_range(*temp_coord, range, abs_min);
            //printf("at [x][y] = [%d][%d] = %f\n", c, r, *temp_coord);

        }
    }
    */
};

//Reading and writing files
int read_bmp(char* image, my_image_comp* input_comps, int* num_comps, int H) {
    //image = argv[1] or arv[2]
    //H is border required
    bmp_in in;
    int err_code = 0;
    if ((err_code = bmp_in__open(&in, image)) != 0) {
        return err_code;
    }

    //extract information on dimension
    int width = in.cols, height = in.rows;
    int n = in.num_components;
    *num_comps = in.num_components;
    input_comps->width = width;
    input_comps->height = height;

    //intiate block of memory
    input_comps->init(height, width, H);

    //copy image pixels into memory block
    int count = 0; //pixel counter
    int r; // Declare row index
    io_byte* line = new io_byte[width * (*num_comps)];
    for (r = height - 1; r >= 0; r--)
    { // "r" holds the true row index we are reading, since the image is
      // stored upside down in the BMP file.
        if ((err_code = bmp_in__get_line(&in, line)) != 0)
            return err_code;
        for (n = 0; n < (*num_comps); n++)
        {
            io_byte* src = line + n; // Points to first sample of component n
            float* dst = input_comps->buf + r * input_comps->stride;
            for (int c = 0; c < width; c++, src += (*num_comps)) {
                dst[c] = (float)*src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
                //count++;
            }
        }
    }
    bmp_in__close(&in);
    //printf("count = %d\n", count);
    delete[] line;
    return 0;
}
float make_in_range(float pattern_val, float range, float abs_min){
    float in_range = ((pattern_val + abs_min)/range)*255.0f;
    return in_range;
}
int write_bmp(my_image_comp* output_comps, char* dest) {
    bmp_out out;
    int num_comps = 1;
    int err_code = 0;
    int height = output_comps->height;
    int width = output_comps->width;


    io_byte* line = new io_byte[width];
    if ((err_code = bmp_out__open(&out, dest, width, height, num_comps)) != 0)
        return err_code;
    for (int r = height - 1; r >= 0; r--)
    { // "r" holds the true row index we are writing, since the image is
      // written upside down in BMP files.
        for (int n = 0; n < num_comps; n++)
        {
            io_byte* dst = line + n; // Points to first sample of component n
            float* src = output_comps->buf + r * output_comps->stride;
            for (int c = 0; c < width; c++, dst += num_comps)
                if (src[c] < 0) {
                    *dst = 0;
                }
                else {
                    *dst = (io_byte)(src[c]); // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
                      //Do I use static_cast? Or do I implement a union?
                }
        }
        bmp_out__put_line(&out, line);
    }
    bmp_out__close(&out);
    delete[] line;
    return 0;
}
/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char* argv[])
{
    if (argc != 4) //exe name, 2 input files, 1 output file
    {
        fprintf(stderr, "Usage: %s <search bmp file> <pattern bmp file> <output bmp file>\n", argv[0]);
        return -1;
    }

    int err_code = 0;
    try {
        //----------------------------------------------------------//
        //                     Design Filter                        //
        //----------------------------------------------------------//
        //4 types of filters, returns centre and length of filter (filter 0 and 1 look ok, the others seem wrong)
        filt filter = make_filter(1);
        //Check filter and store values and dimensions
        int H = (filter.length - 1) / 2;
        int length = filter.length;
        float* hpf = filter.centre;
        for (int r = -H; r <= H; r++) {
            for (int c = -H; c <= H; c++) {
                printf("%f ", *(hpf + r * length + c));
            }
            printf("\n");
        }
        
        //----------------------------------------------------------//
        //                     Read Inputs                          //
        //----------------------------------------------------------//
        my_image_comp search;
        my_image_comp pattern;
        int num_comps = 1;
        //if read_bmp returns -1 ==  fail, then throw err_code which read_bmp will return
        if ((err_code = read_bmp(argv[1], &search, &num_comps, 0)) != 0) {
            throw err_code;
        }
        if ((err_code = read_bmp(argv[2], &pattern, &num_comps, H)) != 0) {
            throw err_code;
        }
        pattern.perform_boundary_extension();
        search.perform_boundary_extension();
        //check values have been extracted properly
        /*
        printf("search:\n");
        printf("    -height = %d\n", search.height);
        printf("    -width  = %d\n", search.width);
        printf("    -stride = %d\n", search.stride);
        printf("pattern:\n");
        printf("    -height = %d\n", pattern.height);
        printf("    -width  = %d\n", pattern.width);
        printf("    -stride = %d\n", pattern.stride);
        */
        //----------------------------------------------------------//
        //                      Processing                          //
        //----------------------------------------------------------//
        //Pattern through HPF: is filter wrong?
        //      -Symmetric
        //make output block
        my_image_comp out;
        out.width = pattern.width;
        out.stride = out.width;
        out.height = pattern.height;
        out.buf = new float[pattern.width * pattern.height];
        /*
        for (int r = 0; r < search.height; r++) {
            for (int c = 0; c < search.width; c++) {
                float* temp_coord = search.buf + r * search.stride + c;
                *(search.buf + r * search.stride + c) = convolve(temp_coord, search.stride, hpf, H);
                if (*temp_coord > max) {
                    max = *temp_coord;
                }
                if (*temp_coord < min) {
                    min = *temp_coord;
                }
            }
        }
        */

        min_max min_max = apply_hpf(pattern.buf, out.buf, pattern.height, pattern.width, pattern.stride, hpf, H);

        
        //Implement Correlation, seems correct
        coord coords = find_coord(&search, &out);
        printf("Correlation at [x][y] = [%d][%d]\n", coords.x, coords.y);
        
        //----------------------------------------------------------//
        //                     Write outputs                        //
        //----------------------------------------------------------//
        //Make filtered pattern between 0 and 255
        float abs_min = fabs(min_max.min);
        float range = min_max.max + abs_min;
        for (int r = 0; r < out.height; r++) {
            for (int c = 0; c < out.width; c++) {
                float* temp_coord = out.buf + r * out.width + c;
                *temp_coord = make_in_range(*temp_coord, range, abs_min);
                //printf("at [x][y] = [%d][%d] = %f\n", c, r, *temp_coord);

            }
        }

        if ((err_code = write_bmp(&out, argv[3])) != 0) {
            throw err_code;
        }

        delete[] out.buf;
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
