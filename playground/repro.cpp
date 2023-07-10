#include <memory>
#include <string>

#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  mbpt::set_default_convention();

  auto dummy = std::make_shared<IndexRegistry>();
  //auto tester = std::make_shared<boost::container::flat_map<int, int>>();
  //auto tester2 = std::make_shared<boost::container::flat_map<int, int>>();

  auto test = std::wstring(L"abcd") < L"bcdef";
}
