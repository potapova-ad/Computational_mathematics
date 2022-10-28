#include "seidel.hpp"

int main()
{
    double da = 0.01;
    double dx = 0.01;
    int n_max = 7;
    int iter_max = 1000;
    std::vector<int> n_line = {};
    std::vector<double> a_line = {};
    std::vector<int> k_line = {};

    for (int n = 1; n <= n_max; ++n)
    {

        for (double a = 0; a <= 1; a += da)
        {
            matrix<double> A = get_matrix(n, a);
            matrix<double> f = get_vector(n, a);
            matrix<double> ans = get_answer(n, a);

            seidel_method seidel(A, f);

            int itter_count = 0;
            double delta = dx * 2;
            do
            {
                if (itter_count == iter_max)
                    break;
                delta = std::abs(norm().n(ans) - norm().n(seidel.get_k_itter(itter_count)));
                ++itter_count;
            } while (delta > dx);
            --itter_count;

            if (itter_count < 1000)
            {
                n_line.push_back(n);
                a_line.push_back(a);
                k_line.push_back(itter_count);
            }
        }
    }

    auto l = matplot::plot3(a_line, n_line, k_line)->line_width(1).color("blue");
    matplot::xlabel("a");
    matplot::ylabel("n");
    matplot::zlabel("Количество необходимых иттераций");
    matplot::show();

    return 0;
}