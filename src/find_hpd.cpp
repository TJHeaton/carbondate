#include <vector>
#include "cpp11.hpp"
using namespace cpp11;

double find_probability_and_ranges_for_cut_off(
        double cut_off,
        doubles calendar_ages,
        doubles density,
        std::vector<double> &start_ages,
        std::vector<double> &end_ages,
        std::vector<double> &auc) {
    start_ages.clear();
    end_ages.clear();
    auc.clear();
    double y1, y2, dx, res = calendar_ages[1] - calendar_ages[0];
    double range_probability = 0, total_probability = 0;
    for (int i = 0; i < calendar_ages.size() - 1; i++) {
        y1 = density[i];
        y2 = density[i + 1];
        if (y1 <= cut_off and y2 > cut_off) {
            dx = res * (cut_off - y1) / (y2 - y1);
            start_ages.push_back(calendar_ages[i] + dx);
            range_probability = (cut_off + y2) * (res - dx) / 2.;
        } else if (y1 > cut_off and y2 <= cut_off) {
            dx = res * (cut_off - y1) / (y2 - y1);
            end_ages.push_back(calendar_ages[i] + dx);
            range_probability += (y1 + cut_off) * dx / 2.;
            auc.push_back(range_probability);
            total_probability += range_probability;
            range_probability = 0;
        } else if (y1 > cut_off and y2 > cut_off) {
            range_probability += (y1 + y2) * res / 2.;
        }
    }
    return total_probability;
}


[[cpp11::register]] list FindHPD(doubles calendar_ages, doubles density, double probability) {

    // Use bisection method to find the closest probability cut-off to give the desired probability
    // a and b are the upper and lower points of the section, p is the midpoint between a and b
    double a = 0., b = 0., p;
    std::vector<double> start_ages, end_ages, auc;
    double current_probability;
    const int max_iter = 1000;

    // Get the maximum value of the density
    for (double dens : density) if (dens > b) b = dens;

    for (int i = 0; i < max_iter; i++) {
        p = (a + b) / 2.;
        current_probability = find_probability_and_ranges_for_cut_off(
            p, calendar_ages, density, start_ages, end_ages, auc);
        if (abs(current_probability - probability) < 1e-4) {
            break;
        }
        if (current_probability < probability) {
            b = p;
        } else {
            a = p;
        }
    }

    writable::list retdata({
    "start_ages"_nm = start_ages,
    "end_ages"_nm = end_ages,
    "area_under_curve"_nm = auc,
    "height"_nm = p,
    });

    return retdata;
}