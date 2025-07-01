#ifndef SEQUANT_CORE_EXPORT_TAPP_HPP
#define SEQUANT_CORE_EXPORT_TAPP_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/reordering_context.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <range/v3/view/enumerate.hpp>

#include <cassert>
#include <format>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>

namespace sequant {

class TappContext : public ExportContext {
 public:
  TappContext() = default;
};

template <typename Context = TappContext>
class TappGenerator : public Generator<Context> {
 public:
  TappGenerator() = default;
  ~TappGenerator() = default;

  std::string get_format_name() const override { return "TAPP"; }

  bool supports_named_sections() const override { return false; }

  bool requires_named_sections() const override { return false; }

  DeclarationScope index_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope variable_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  DeclarationScope tensor_declaration_scope() const override {
    return DeclarationScope::Global;
  }

  std::string represent(const Index &idx, const Context &ctx) const override {
	  assert(toUtf8(idx.label()).size() == 1);
    return toUtf8(idx.label());
  }

  std::string represent(const Tensor &tensor,
                        const Context &ctx) const override {
    std::string representation = toUtf8(label(tensor));

	if (tensor.num_indices() > 0) {
		representation += "_";
		for (const Index &idx : tensor.indices()) {
			representation += toUtf8(idx.space().base_key());
		}
	}

    return representation;
  }

  std::string represent(const Variable &variable,
                        const Context &ctx) const override {
    return toUtf8(variable.label());
  }

  std::string represent(const Constant &constant,
                        const Context &ctx) const override {
    std::stringstream sstream;
    if (constant.value().imag() != 0) {
		throw std::string("Complex-valued constants not (yet) supported");
    } else {
      sstream << constant.value().real();
    }

    return sstream.str();
  }

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx) override {
    if (!zero_init) {
		// TODO
    }
  }

  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx) override {
  }

  void set_to_zero(const Tensor &tensor, const Context &ctx) override {
  }

  void unload(const Tensor &tensor, const Context &ctx) override {
  }

  void destroy(const Tensor &tensor, const Context &ctx) override {
  }

  void persist(const Tensor &tensor, const Context &ctx) override {
  }

  void create(const Variable &variable, bool zero_init,
              const Context &ctx) override {
	  m_generated += "float " + represent(variable, ctx);
  }

  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx) override {
  }

  void set_to_zero(const Variable &variable, const Context &ctx) override {
  }

  void unload(const Variable &variable, const Context &ctx) override {
  }

  void destroy(const Variable &variable, const Context &ctx) override {
  }

  void persist(const Variable &variable, const Context &ctx) override {
  }
  

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx) override {
	  // TODO
  }

  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx) override {
	  // TODO
  }

  void declare(const Index &idx, const Context &ctx) override {
  }

  void declare(const Variable &variable, UsageSet usage,
               const Context &ctx) override {
  }

  void declare(const Tensor &tensor, UsageSet usage,
               const Context &ctx) override {
	std::string extents = "{ ";
	std::string strides = "{ ";
	std::string next_stride = "1";
	const std::size_t num_indices = tensor.num_indices();
	for (const auto [i, idx] : ranges::views::enumerate(tensor.indices())) {
		std::string current = "num_" + toUtf8(idx.space().base_key());
		extents += current;
		strides += next_stride;

		if (i == 0) {
			next_stride = current;
		} else {
			next_stride += " * " + current;
		}

		if (i + 1 < num_indices) {
			extents += ", ";
			strides += ", ";
		}
	}

	extents += "}";
	strides += "}";

	std::string info = std::format("info_{}", represent(tensor, ctx));
	m_generated += std::format("TAPP_tensor_info {};\n", info);
	m_generated += std::format("TAPP_create_tensor_info(&{}, TAPP_F32, {}, {}, {});\n", info, num_indices, extents, strides);
  }

  void all_indices_declared(std::size_t amount, const Context &ctx) override {
  }

  void all_variables_declared(std::size_t amount, const Context &ctx) override {
  }

  void all_tensors_declared(std::size_t amount, const Context &ctx) override {
  }

  void begin_declarations(DeclarationScope scope, const Context &ctx) override {
  }

  void end_declarations(DeclarationScope scope, const Context &ctx) override {
  }

  void insert_comment(const std::string &comment, const Context &ctx) override {
    m_generated += "// " + comment + "\n";
  }

  void begin_named_section(std::string_view name, const Context &ctx) override {
  }

  void end_named_section(std::string_view name, const Context &ctx) override {
  }

  void begin_expression(const Context &ctx) override {
  }

  void end_expression(const Context &ctx) override {}

  void begin_export(const Context &ctx) override { m_generated.clear(); }

  void end_export(const Context &ctx) override { }

  std::string get_generated_code() const override { return m_generated; }

 private:
  std::string m_generated;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_TAPP_HPP
