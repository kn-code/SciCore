//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#ifndef SCICORE_SERIALIZATION_H
#define SCICORE_SERIALIZATION_H

#include "Definitions.h"
#include "SciCore_export.h"

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>

namespace SciCore::Detail
{

// We don't want to collide with the default cereal implementation
// for complex numbers. We therefore make our own wrapped complex
// number type.
// The reason not to use the default implementation is that it saves
// every complex number with a named "real" and "imag" field, which
// leads to a large size overhead for complex matrices.
struct SCICORE_EXPORT WrappedComplex
{
    Complex c;
};
} // namespace SciCore::Detail

namespace cereal
{

template <typename Archive, traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void save(Archive& archive, const SciCore::Detail::WrappedComplex& x)
{
    cereal::size_type s = 2;
    archive(cereal::make_size_tag(s));
    archive(x.c.real(), x.c.imag());
}

template <typename Archive, traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void load(Archive& archive, SciCore::Detail::WrappedComplex& x)
{
    cereal::size_type s;
    archive(cereal::make_size_tag(s));
    if (s != 2)
    {
        throw std::runtime_error("load WrappedComplex failed: wrong format");
    }
    SciCore::Real real, imag;
    archive(real, imag);
    x.c.real(real);
    x.c.imag(imag);
}

// std::vector<Complex> -- text archive
template <typename Archive, traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void save(Archive& archive, const std::vector<SciCore::Complex>& x)
{
    cereal::size_type s = x.size();
    archive(cereal::make_size_tag(s));

    for (const auto& y : x)
    {
        SciCore::Detail::WrappedComplex w{y};
        archive(w);
    }
}

template <typename Archive, traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void load(Archive& archive, std::vector<SciCore::Complex>& x)
{
    cereal::size_type s;
    archive(cereal::make_size_tag(s));

    x.resize(s);

    for (size_t i = 0; i < x.size(); ++i)
    {
        SciCore::Detail::WrappedComplex w;
        archive(w);
        x[i] = w.c;
    }
}

// std::vector<Complex> -- binary archive
template <typename Archive, traits::EnableIf<traits::is_output_serializable<BinaryData<SciCore::Complex>, Archive>::value> = traits::sfinae>
SCICORE_EXPORT void save(Archive& archive, const std::vector<SciCore::Complex>& x)
{
    size_t size = x.size();
    archive(size);
    archive(binary_data(x.data(), size * sizeof(SciCore::Complex)));
}

template <typename Archive, traits::EnableIf<traits::is_input_serializable<BinaryData<SciCore::Complex>, Archive>::value> = traits::sfinae>
SCICORE_EXPORT void load(Archive& archive, std::vector<SciCore::Complex>& x)
{
    size_t size;
    archive(size);

    x.resize(size);

    archive(binary_data(x.data(), size * sizeof(SciCore::Complex)));
}

// Column vector
template <
    typename Archive,
    class _Scalar,
    int _Rows,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void save(Archive& archive, const Eigen::Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols>& m)
{
    cereal::size_type s = m.rows();
    archive(cereal::make_size_tag(s));

    for (const auto& x : m)
    {
        if constexpr (std::is_same_v<SciCore::Real, std::decay_t<_Scalar>> == true)
        {
            archive(x);
        }
        else // Complex
        {
            SciCore::Detail::WrappedComplex wx{x};
            archive(wx);
        }
    }
}

template <
    typename Archive,
    class _Scalar,
    int _Rows,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void load(Archive& archive, Eigen::Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols>& m)
{
    cereal::size_type s;
    archive(cereal::make_size_tag(s));

    m.resize(s);

    for (int i = 0; i < m.size(); ++i)
    {
        if constexpr (std::is_same_v<SciCore::Real, std::decay_t<_Scalar>> == true)
        {
            archive(m[i]);
        }
        else // Complex
        {
            SciCore::Detail::WrappedComplex wx;
            archive(wx);
            m[i] = wx.c;
        }
    }
}

// Row vector
template <
    typename Archive,
    class _Scalar,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void save(Archive& archive, const Eigen::Matrix<_Scalar, 1, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
    cereal::size_type s = m.cols();
    archive(cereal::make_size_tag(s));

    for (const auto& x : m)
    {
        if constexpr (std::is_same_v<SciCore::Real, std::decay_t<_Scalar>> == true)
        {
            archive(x);
        }
        else // Complex
        {
            SciCore::Detail::WrappedComplex wx{x};
            archive(wx);
        }
    }
}

template <
    typename Archive,
    class _Scalar,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
SCICORE_EXPORT void load(Archive& archive, Eigen::Matrix<_Scalar, 1, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
    cereal::size_type s;
    archive(cereal::make_size_tag(s));

    m.resize(s);

    for (int i = 0; i < m.size(); ++i)
    {
        if constexpr (std::is_same_v<SciCore::Real, std::decay_t<_Scalar>> == true)
        {
            archive(m[i]);
        }
        else // Complex
        {
            SciCore::Detail::WrappedComplex wx;
            archive(wx);
            m[i] = wx.c;
        }
    }
}

// Matrices
template <
    typename Archive,
    class _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
    requires(_Rows != 1 || _Cols != 1)
SCICORE_EXPORT void save(Archive& archive, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
    int rows = m.rows();

    archive(cereal::make_size_tag(rows));
    for (int i = 0; i < rows; ++i)
    {
        Eigen::Matrix<_Scalar, -1, 1> row = m.row(i);
        archive(row);
    }
}

template <
    typename Archive,
    class _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_text_archive<Archive>::value> = traits::sfinae>
    requires(_Rows != 1 || _Cols != 1)
SCICORE_EXPORT void load(Archive& archive, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
    size_t rows;
    archive(cereal::make_size_tag(rows));

    using DynamicVec = Eigen::Matrix<_Scalar, -1, 1>;
    std::vector<DynamicVec> rowEntries(rows);
    for (size_t i = 0; i < rows; ++i)
    {
        archive(rowEntries[i]);
    }

    if (rows == 0)
    {
        m = Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>();
        return;
    }

    int cols = rowEntries[0].size();
    m.resize((int)rows, cols);
    for (int i = 0; i < (int)rows; ++i)
    {
        if (rowEntries[i].size() != cols)
        {
            throw std::runtime_error("load Matrix failed: wrong number of columns");
        }

        m.row(i) = rowEntries[i];
    }
}

// Save & load into/from binary archives
template <
    class Archive,
    class _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value> = traits::sfinae>
SCICORE_EXPORT void save(Archive& archive, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
    int rows = m.rows();
    int cols = m.cols();
    archive(rows);
    archive(cols);
    archive(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
}

template <
    class Archive,
    class _Scalar,
    int _Rows,
    int _Cols,
    int _Options,
    int _MaxRows,
    int _MaxCols,
    traits::EnableIf<traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value> = traits::sfinae>
SCICORE_EXPORT void load(Archive& archive, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
    int rows;
    int cols;
    archive(rows);
    archive(cols);

    m.resize(rows, cols);

    archive(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
}

} // namespace cereal

#endif // SCICORE_SERIALIZATION_H
