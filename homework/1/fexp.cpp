#include "maclaurin.hpp"
using namespace matplot;

Maclaurin::Maclaurin(const function_type &f)
{
    _der = new forward_der(f);

}


int main()
{
    auto f = [](double x) { return std::exp(x); };

    Maclaurin f0(f);

    double xbegin = 0, xend = 1, dt = 0.001;
    int kbegin = 1, dk = 1;

    int akmax = 6;


        std::vector<int> x = {};
        std::vector<double> y = {};

        for (int k = kbegin; k <= akmax; k += dk)
        {
            double sum = 0;
            int count = 0;

            for (double t = xbegin; t <= xend; t += dt)
            {
                sum += std::abs(f0.value(k, t, dt) - f(t));
                ++count;
            }
            x.push_back(k);
            y.push_back(sum / count);
        }

        title("Производная вперед для exp(x)");
        grid(true);
        xlabel("n");
        ylabel("Среднее отклонение");
        plot(x, y)->line_width(1).color("blue");
        show();
    

    return 0;
}