#include <fstream>
#include "transformer.h"

int main()
{
    std::ifstream stream("../data/highway_map.csv");
    transformer_t transformer(stream);

    auto x = transformer(6914.15);
    std::cout << x.x << " " << x.y << "\n";

    return 0;
}
