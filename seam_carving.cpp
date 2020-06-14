#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class Image {
public:    

    int W;
    int H;
    vector<int> colors; // W * H * 3
    vector<int> beam; // all x values forming the vertical beam
    
    Image(int width, int height) {
        W = width;
        H = height;
        colors.resize(W * H);
    }

    int intensity(int x, int y) {

        return colors[(x * W + y)*3 + 0] 
             + colors[(x * W + y)*3 + 1]
             + colors[(x * W + y)*3 + 2];

    }

    int energyMap(int x, int y) {
        
        int up_intensity = (y == 0) ? 0 : intensity(x, y + 1);
        int left_intensity = (x == 0) ? 0 : intensity(x - 1, y);
        int right_intensity = (x == W - 1) ? 0 : intensity(x + 1, y);
        int down_intensity = (y == H - 1) ? 0 : intensity(x, y - 1);

        return abs(left_intensity - right_intensity) + abs(up_intensity - down_intensity);

    }

    int valueFunction(int x, int y) {

        int cur_energy = energyMap(x, y);
    
        if (y == 0)
            return cur_energy;

        std::vector<int> nxt_nrg = {valueFunction(x - 1, y - 1), valueFunction(x, y - 1), valueFunction(x + 1, y - 1)};

        auto it = min_element(nxt_nrg.begin(), nxt_nrg.end());
        int min_idx = it - nxt_nrg.begin();

        if (min_idx == 0)
            beam.push_back(x - 1);
        
        else if (min_idx == 1)
            beam.push_back(x);

        else
            beam.push_back(x + 1);
        
        return *it + cur_energy;

    }

};

int main() {

    int W = 512, H = 512;
    Image img(W, H);

    int pixels_to_remove = 50;

    for (int t = 0; t < pixels_to_remove; t++){

        // find lowest energy cell in last row
        vector<int> energies;
        for (int x = 0; x < img.W; ++x) {
            int cell_energy = img.energyMap(x, img.H - 1);
            energies.push_back(cell_energy);
        }

        int min_x = min_element(energies.begin(), energies.end()) - energies.begin();

        // now backtrace from (min_x, H - 1)
        img.beam.push_back(min_x);
        img.valueFunction(min_x, img.H - 1);

        // remove beam from img (i.e. copy non-beam pixels to new img)
        Image new_img(img.W - 1, img.H);
        for (int i = 0; i < img.W * img.H * 3; ++i) {

            int y = (int) i / img.W;
            if (i % img.W == img.beam[img.H - y - 1])
                continue;

            new_img.colors.push_back(img.colors[i]);
        }

        img = new_img;

    }

    return 0;

}