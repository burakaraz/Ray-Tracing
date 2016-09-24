#include <cstdio>

void writePPM(const char* filename, int width, int height, float *data)
{
    FILE* fp = fopen(filename, "w");

    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");

    for (int y = 0, i = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x, i += 3)
        {
            fprintf(fp, "%d %d %d ", (int) data[i], (int) data[i+1], (int) data[i+2]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
