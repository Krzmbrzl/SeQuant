#include <peglib.h>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/attr.hpp>

#include <optional>
#include <string_view>

namespace sequant {

ResultExpr parse_result_expr(std::wstring_view input,
                             std::optional<Symmetry> perm_symm,
                             std::optional<BraKetSymmetry> braket_symm,
                             std::optional<ColumnSymmetry> column_symm) {
}

}
