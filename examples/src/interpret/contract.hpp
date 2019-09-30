#ifndef SEQUANT_INTERPRET_CONTRACT_HPP
#define SEQUANT_INTERPRET_CONTRACT_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <btas/btas.h>
#include <btas/tensorview.h>
#include <btas/tensor_func.h>

#include <tiledarray.h>

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>

#include "interpreted_tensor.hpp"

namespace sequant {

  namespace interpret {

    template <typename T>
      InterpretedTensor<T> eval_equation(const sequant::ExprPtr&,
          const std::map<std::wstring, T const * >&);

    template <typename T>
      void antisymmetrize(InterpretedTensor<T>&, const size_t&);

    namespace detail {

      template <typename T>
        InterpretedTensor<T> contract(const InterpretedTensor<T>&,
            const InterpretedTensor<T>&, double scal=1.0);

      template <typename T>
        InterpretedTensor<T> contract_vec(const std::vector<InterpretedTensor<T>>&,
            size_t, double scal=1.0, size_t b=1);

      template <typename T>
        InterpretedTensor<T> eval_product(const sequant::ExprPtr&,
            const std::map<std::wstring, T const * >&);

      template <typename T>
        InterpretedTensor<T> eval_sum(const sequant::ExprPtr&,
            const std::map<std::wstring, T const * >&);

      template <typename T>
        InterpretedTensor<T> sum(const InterpretedTensor<T>&,
            const InterpretedTensor<T>&, bool subtract=false);

      const char ADD{ 1}, SUB{-1};
      struct Perm { std::vector<size_t> perm; char even_perm{ADD}; };

      std::vector<Perm> perm_calc(std::vector<size_t>,
          size_t,
          size_t cswap = 0, // count swaps
          size_t begin = 0);

      auto permute(const btas::Tensor<double>&, const std::vector<size_t>&);

      auto permute(const TA::TArrayD&, const std::vector<size_t>&);

      std::string ords_to_csv_str(const std::vector<size_t>&);

      std::string range_to_csv_str(const size_t&);

      btas::Tensor<double> core_contract(double, // scalar
          const btas::Tensor<double>&, // first tensor
          const std::vector<size_t>&,  // indexes of the first tensor
          const btas::Tensor<double>&, // second tensor
          const std::vector<size_t>&,  // indexes of the second tensor
          const std::vector<size_t>&   // non-contracting indices
          );

      TA::TArrayD core_contract(double,
          const TA::TArrayD&,
          const std::vector<size_t>&,
          const TA::TArrayD&,
          const std::vector<size_t>&,
          const std::vector<size_t>&
          );

      btas::Tensor<double> scale(double, const btas::Tensor<double>&);

      TA::TArrayD scale(double, const TA::TArrayD&);

      btas::Tensor<double> core_sum(const btas::Tensor<double>&,
          const btas::Tensor<double>&, bool subtract=false);

      TA::TArrayD core_sum(const TA::TArrayD&,
          const TA::TArrayD&, bool subtract=false);

    } // namespace detail

  } // namespace interpret

} // namespace sequant

namespace sequant {

  namespace interpret {

    template <typename T>
    InterpretedTensor<T> eval_equation(const sequant::ExprPtr& expr,
        const std::map<std::wstring, T const * >& tmap) {

      if (expr->is<sequant::Product>())
        return detail::eval_product(expr, tmap);

      else if (expr->is<sequant::Sum>())
        return detail::eval_sum(expr, tmap);

      else if (expr->is<sequant::Tensor>()) {
        auto ct = InterpretedTensor<T>(expr->as<sequant::Tensor>());
        ct.link_tensor(*(tmap.find(ct.label()+L"_"+ct.translate())->second));
        return ct;
      }

      else
        throw "Only know how to handle sum or product!";
    }

    template <typename T>
    void antisymmetrize(InterpretedTensor<T>& tensor, const size_t& rank) {
      if (rank == 2)
        return;
      else if (rank%2 != 0)
        throw "Can't handle odd-ordered tensors, just yet!";

      std::vector<size_t> to_perm;
      for (auto i = 0; i < (size_t)rank/2; ++i)
        to_perm.push_back(i);
      auto vp = detail::perm_calc(to_perm, (size_t)rank/2);

      // antisymmetrize
      auto result = tensor.tensor();
      result.fill(0); // could prove fragile
                      // coincidentally both btas::Tensor and TA::DistArray
                      // objects have a fill method with the same name

      for (const auto& p: vp) {
        for (const auto& q: vp) {
          // permutation of the bra
          auto perm_vec = p.perm;
          // permutation of the ket
          // ket indices = rank/2 + bra indices
          for (auto qq: q.perm)
            perm_vec.push_back((size_t)rank/2  + qq);

          auto perm_t = detail::permute(tensor.tensor(), perm_vec);

          if (p.even_perm * q.even_perm == detail::ADD)
            result = detail::core_sum(result, perm_t);
          else
            result = detail::core_sum(result, perm_t, true); // subtract = true
        } // for q: vp
      } // for p: vp
      tensor.link_tensor(result);
    } // function antisymmetrize

    namespace detail {

      template <typename T>
        InterpretedTensor<T> contract(const InterpretedTensor<T>& t1,
            const InterpretedTensor<T>& t2, double scal) {

          using ords_vec    = std::vector<size_t>;
          using index_vec   = std::vector<sequant::Index>;

          using ind_ord     = std::pair<sequant::Index, size_t>;
          using ind_ord_map = std::map<ind_ord::first_type, ind_ord::second_type>;

          ind_ord_map index_to_ordinal; // all indices to ordinal map
          ind_ord_map nc_indices;       // non-contracting indices to ordinal map
          size_t count = 0;

          auto update_index_to_ordinal = [&](const index_vec& vec){

            auto map_size = index_to_ordinal.size();

            for (const auto& i: vec) {
              index_to_ordinal.insert(ind_ord(i, ++count));
              auto current_size = index_to_ordinal.size();
              if (map_size < current_size) {
                nc_indices.insert(ind_ord(i, count));
                ++map_size;
              }
              else
                nc_indices.erase(i);
            }
          };

          update_index_to_ordinal(t1.bras());
          update_index_to_ordinal(t1.kets());
          update_index_to_ordinal(t2.bras());
          update_index_to_ordinal(t2.kets());

          auto gen_ords_vec = [&index_to_ordinal](const InterpretedTensor<T>& t){
            ords_vec v;
            for (const auto& i: t.bras())
              v.push_back(index_to_ordinal.at(i));

            for (const auto& i: t.kets())
              v.push_back(index_to_ordinal.at(i));
            return v;
          };

          const ords_vec t1_ords = gen_ords_vec(t1);

          const ords_vec t2_ords = gen_ords_vec(t2);

          ords_vec nc_ords;
          index_vec nc_index_vec;
          for (const auto& i: nc_indices) {
            nc_index_vec.push_back(i.first);
            nc_ords.push_back(i.second);
          }

          InterpretedTensor<T> result{t1.label(), nc_index_vec};

          result.link_tensor(core_contract(scal, t1.tensor(),
                t1_ords, t2.tensor(), t2_ords, nc_ords));

          return result;
        }

      template <typename T>
        InterpretedTensor<T> contract_vec(const std::vector<InterpretedTensor<T>>& vct,
            size_t n, double scal, size_t b) {
          // evaluates left to right
          auto i = n - b;
          if (i == 0) {
            if (scal != 1.0) {
              auto tmp = vct[i];
              auto tmp_tnsr = vct[i].tensor();
              tmp.link_tensor(scale(scal, tmp_tnsr));
              return tmp;
            }
            return vct[i];
          }
          else if (i == 1)
            return contract(vct[i-1], vct[i], scal);
          // else
          return contract(vct[n-b], contract_vec(vct, n, 1.0, b+1), scal);
        }

      template <typename T>
        InterpretedTensor<T> eval_product(const sequant::ExprPtr& expr,
            const std::map<std::wstring, T const * >& tmap) {
          auto p = expr->as<sequant::Product>();

          // collect the tensors to be contracted
          std::vector<InterpretedTensor<T>> fvec;
          for (const auto& f: p.factors()) {
            if (!f->is<sequant::Tensor>())
              // call equation evaluater if factor is not a tensor
              fvec.push_back(eval_equation(f, tmap));
            else {
              auto t = f->as<sequant::Tensor>();
              // "A": antisymmetrizing tensors are ommitted at this point
              if (t.label() != L"A")
                fvec.push_back(InterpretedTensor<T>{t});
            }
          }
          // linking the suitable tensors from the map
          for (auto& f: fvec) {
            f.link_tensor(*(tmap.find(f.label() + L"_" + f.translate())->second));
          }
          return contract_vec(fvec, fvec.size(), std::real(p.scalar()));
        }

      template <typename T>
        InterpretedTensor<T> eval_sum(const sequant::ExprPtr& expr,
            const std::map<std::wstring, T const * >& tmap) {

          auto smands = expr->as<sequant::Sum>().summands();
          // collect the tensors to be summed
          std::vector<InterpretedTensor<T>> svec;
          //
          for (auto i=0; i < smands.size(); ++i) {
            if (! smands.at(i)->is<sequant::Tensor>()) 
              // call equation evaluater if summand is not a sequant::tensor
              svec.push_back(eval_equation(smands[i], tmap));
            else {
              auto t = smands[i]->as<sequant::Tensor>();
              if (t.label() != L"A") { // 'A' tensors are omitted
                svec.push_back(InterpretedTensor<T>{t});
                svec[i].link_tensor(*(tmap.find(svec[i].label() +
                        L"_" + svec[i].translate())->second));
              }
            }
          }
          auto result = svec[0];
          for (auto i = 1; i < svec.size(); ++i) {
            result = sum(result, svec[i]);
          }
          return result;
        }

      template <typename T>
        InterpretedTensor<T> sum(const InterpretedTensor<T>& s1,
            const InterpretedTensor<T>& s2, bool subtract) {
          auto result = s1; // or s2
          auto t = core_sum(s1.tensor(), s2.tensor(), subtract);
          result.link_tensor(t);
          return result;
        }

      std::vector<Perm> perm_calc(std::vector<size_t> to_perm,
          size_t  size, size_t cswap, size_t begin) {

        if (begin+1 == size) {
          // even permutations added. odds subtracted
          return std::vector<Perm>{Perm{to_perm, (cswap%2==0)? ADD : SUB}};
        }

        auto result = std::vector<Perm>{};
        for (auto i = begin; i < size; ++i) {
          std::swap(to_perm[begin], to_perm[i]);
          auto more_result = perm_calc(to_perm, size,
              (begin==i)? cswap: cswap+1, begin+1);
          for (auto p: more_result) {
            result.push_back(p);
          }
          std::swap(to_perm[begin], to_perm[i]);
        }
        return result;
      }

      auto permute(const btas::Tensor<double>& bt,
          const std::vector<size_t>& perm_vec) {
        return btas::Tensor<double>(btas::permute(bt, perm_vec));
      }


      auto permute(const TA::TArrayD& ta,
          const std::vector<size_t>& perm_vec) {

        // assumes we always get the perm_vec whose elements
        // are the contiguous elements from 0 to n-1 where
        // n is the size of the perm_vec

        auto permed = ta;
        permed.fill(0);
        permed(range_to_csv_str(perm_vec.size())) = permed(ords_to_csv_str(perm_vec));
        return permed;
      }

      std::string ords_to_csv_str(const std::vector<size_t>& ords) {

        std::string str = "";
        for (const auto& o: ords)
          str += std::to_string(o) + ",";

        str.pop_back(); // remove the trailing comma(,)
        return str;
      }

      std::string range_to_csv_str(const size_t& n) {
        std::vector<size_t> range_vec(n);
        for (auto i = 0; i < n; ++i)
          range_vec[i] = i;
        return ords_to_csv_str(range_vec);
      }

      btas::Tensor<double> core_contract(double scal, // scalar
          const btas::Tensor<double>& t1,      // first tensor
          const std::vector<size_t>&  t1_ords, // indexes of the first tensor
          const btas::Tensor<double>& t2,      // second tensor
          const std::vector<size_t>&  t2_ords, // indexes of the second tensor
          const std::vector<size_t>&  nc_ords  // non-contracting indices
          ) {
        btas::Tensor<double> result;
        btas::contract(scal,
            t1, t1_ords, t2, t2_ords,
            0.0, result, nc_ords);
        return result;
      }

      TA::TArrayD core_contract(double scal,
          const TA::TArrayD& t1,
          const std::vector<size_t>& t1_ords,
          const TA::TArrayD& t2,
          const std::vector<size_t>& t2_ords,
          const std::vector<size_t>& nc_ords) {
        TA::TArrayD result;

        result(ords_to_csv_str(nc_ords))
          = t1(ords_to_csv_str(t1_ords))
          * t2(ords_to_csv_str(t2_ords));
        return result;
      }

      btas::Tensor<double> scale(double d, const btas::Tensor<double>& t) {
        auto result = t;
        btas::scal(d, result);
        return result;
      }

      TA::TArrayD scale(double d, const TA::TArrayD& t) {
        auto result = t;
        TA::scale(result, d);
        return result;
      }

      btas::Tensor<double> core_sum(const btas::Tensor<double>& t1,
          const btas::Tensor<double>& t2, bool subtract) {
        if (subtract)
          return t1 - t2;
        // else
        return t1 + t2;
      }

      TA::TArrayD core_sum(const TA::TArrayD& t1,
          const TA::TArrayD& t2, bool subtract) {
        if (!t2.is_initialized())
          return t1;
        if (!t1.is_initialized()) {
          if (!subtract)
            return t2;
          else {
            TArrayD result;
            auto rank    = t1.trange().rank();
            auto inds    = range_to_csv_str(rank);
            result(inds) = t2(inds) * -1;
            return result;
          }
        }
        // now both t1 and t2 are initialized
        auto rank    = t1.trange().rank();
        auto inds    = range_to_csv_str(rank);
        TArrayD result;
        if (subtract)
          result(inds) = t1(inds) - t2(inds);
        else
          result(inds) = t1(inds) + t2(inds);
        return result;
      }

      } // namespace detail


  } // namespace interpret

} // namespace sequant

#endif /* ifndef SEQUANT_INTERPRET_CONTRACT_HPP */
