//
// Created by Eduard Valeyev on 2019-02-18.
//

#include <SeQuant/core/op.hpp>

namespace sequant {
namespace detail {
OpIdRegistrar::OpIdRegistrar() {
  auto id = std::numeric_limits<Expr::type_id_type>::max();
  Expr::set_type_id<FNOperator>(id);
  Expr::set_type_id<BNOperator>(--id);
  Expr::set_type_id<FOperator>(--id);
  Expr::set_type_id<BOperator>(--id);
}
}  // namespace detail

template <>
const container::vector<std::wstring>&
NormalOperator<Statistics::FermiDirac>::labels() {
  using namespace std::literals;
  static container::vector<std::wstring> labels_{L"a"s, L"ã"s};
  return labels_;
}

template <>
const container::vector<std::wstring>&
NormalOperator<Statistics::BoseEinstein>::labels() {
  using namespace std::literals;
  static container::vector<std::wstring> labels_{L"b"s, L"ᵬ"s};
  return labels_;
}

template <>
std::wstring_view NormalOperator<Statistics::FermiDirac>::label() const {
  static const std::wstring a = L"a";
  static const std::wstring atilde = L"ã";
  return vacuum() == Vacuum::Physical ? a : atilde;
}

template <>
std::wstring_view NormalOperator<Statistics::BoseEinstein>::label() const {
  static const std::wstring b = L"b";
  static const std::wstring btilde = L"ᵬ";
  return vacuum() == Vacuum::Physical ? b : btilde;
}

}  // namespace sequant
