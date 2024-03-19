//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

///
/// \file   DCT.h
///
/// \brief  Fast discrete cosine transform.
///

#ifndef SCICORE_DCT_H
#define SCICORE_DCT_H

#include "Definitions.h"

namespace SciCore
{

///
///  \defgroup DCT Fast discrete cosine transform
///
///  \brief One and two dimensional DCTs of scalar and matrix valued data.
///
///  This module contains methods related to DCTs of scalar and matrix valued data.
///
///  \{
///

#define FMA(a, b, c) (((a) * (b)) + (c))
#define FMS(a, b, c) (((a) * (b)) - (c))
#define FNMA(a, b, c) (-(((a) * (b)) + (c)))
#define FNMS(a, b, c) ((c) - ((a) * (b)))

// data and buffer must have 3 elements
template <MatrixOrScalarType T>
void dct_3(T* data, T* buffer, int stride, int bufferStride)
{
    buffer[bufferStride * 0] = data[stride * 1];
    buffer[bufferStride * 1] = data[stride * 2];
    buffer[bufferStride * 2] = buffer[bufferStride * 1] + data[stride * 0];
    data[stride * 1] =
        Real(0.866025403784438646763723170752936183471402627L) * (data[stride * 0] - buffer[bufferStride * 1]);
    data[stride * 0] = buffer[bufferStride * 0] + buffer[bufferStride * 2];
    data[stride * 2] = Real(0.5) * FNMS(Real(2), buffer[bufferStride * 0], buffer[bufferStride * 2]);
}

// data and buffer must have 6 elements
template <MatrixOrScalarType T>
void dct_6(T* data, T* buffer, int stride, int bufferStride)
{
    buffer[bufferStride * 0] = data[stride * 5] + data[stride * 2];
    buffer[bufferStride * 1] = data[stride * 5] - data[stride * 2];
    buffer[bufferStride * 2] = FNMS(Real(0.5), buffer[bufferStride * 0], data[stride * 1]);
    buffer[bufferStride * 3] = data[stride * 4];
    buffer[bufferStride * 4] = data[stride * 0] + data[stride * 3];
    buffer[bufferStride * 5] = data[stride * 0] - data[stride * 3];
    data[stride * 5]         = FMS(Real(0.5), buffer[bufferStride * 4], buffer[bufferStride * 3]);

    data[stride * 2] =
        Real(0.866025403784438646763723170752936183471402627L) * (buffer[bufferStride * 1] + buffer[bufferStride * 5]);
    data[stride * 4] = data[stride * 5] - buffer[bufferStride * 2];

    buffer[bufferStride * 0] += data[stride * 1];
    data[stride * 1]          = buffer[bufferStride * 3] + buffer[bufferStride * 4];
    data[stride * 0]          = buffer[bufferStride * 0] + data[stride * 1];
    data[stride * 3] =
        Real(0.707106781186547524400844362104849039284835938L) * (data[stride * 1] - buffer[bufferStride * 0]);
    buffer[bufferStride * 0] =
        Real(0.612372435695794524549321018676472847991486870L) * (buffer[bufferStride * 5] - buffer[bufferStride * 1]);
    buffer[bufferStride * 3] =
        Real(0.707106781186547524400844362104849039284835938L) * (buffer[bufferStride * 2] + data[stride * 5]);
    data[stride * 5] = buffer[bufferStride * 0] - buffer[bufferStride * 3];
    data[stride * 1] = buffer[bufferStride * 0] + buffer[bufferStride * 3];
}

// data and buffer must have 8 elements
template <MatrixOrScalarType T>
void dct_8(T* data, T* buffer, int stride, int bufferStride)
{
    buffer[bufferStride * 0] = data[stride * 0] + data[stride * 7];
    buffer[bufferStride * 1] = data[stride * 1] + data[stride * 6];
    buffer[bufferStride * 2] = data[stride * 2] + data[stride * 5];
    buffer[bufferStride * 3] = data[stride * 3] + data[stride * 4];
    buffer[bufferStride * 4] = data[stride * 3] - data[stride * 4];
    buffer[bufferStride * 5] = data[stride * 2] - data[stride * 5];
    buffer[bufferStride * 6] = data[stride * 1] - data[stride * 6];
    buffer[bufferStride * 7] = data[stride * 0] - data[stride * 7];

    data[stride * 1]         = buffer[bufferStride * 1] - buffer[bufferStride * 2];
    data[stride * 3]         = buffer[bufferStride * 0] + buffer[bufferStride * 3];
    data[stride * 7]         = buffer[bufferStride * 1] + buffer[bufferStride * 2];
    buffer[bufferStride * 1] = buffer[bufferStride * 0] - buffer[bufferStride * 3];
    buffer[bufferStride * 0] = -buffer[bufferStride * 4] - buffer[bufferStride * 5];
    buffer[bufferStride * 2] = Real(0.707106781186547524400844L) * (buffer[bufferStride * 5] + buffer[bufferStride * 6]);
    buffer[bufferStride * 3] = buffer[bufferStride * 6] + buffer[bufferStride * 7];

    data[stride * 0]         = data[stride * 3] + data[stride * 7];
    data[stride * 4]         = Real(0.707106781186547524400844L) * (data[stride * 3] - data[stride * 7]);
    buffer[bufferStride * 6] = Real(0.707106781186547524400844L) * (data[stride * 1] + buffer[bufferStride * 1]);
    data[stride * 3]         = Real(0.382683432365089771728460L) * (buffer[bufferStride * 0] + buffer[bufferStride * 3]);

    data[stride * 7]         = -Real(0.541196100146196984399723L) * buffer[bufferStride * 0] - data[stride * 3];
    buffer[bufferStride * 0] = Real(1.306562964876376527856643L) * buffer[bufferStride * 3] - data[stride * 3];

    data[stride * 2]         = Real(0.541196100146196984399724L) * (buffer[bufferStride * 6] + buffer[bufferStride * 1]);
    data[stride * 6]         = Real(1.306562964876376527856644L) * (buffer[bufferStride * 1] - buffer[bufferStride * 6]);
    buffer[bufferStride * 1] = buffer[bufferStride * 2] + buffer[bufferStride * 7];
    data[stride * 1]         = buffer[bufferStride * 7] - buffer[bufferStride * 2];

    data[stride * 3] = Real(0.601344886935045280543722L) * (data[stride * 1] - data[stride * 7]);
    data[stride * 5] = Real(0.89997622313641570463851L) * (data[stride * 7] + data[stride * 1]);
    data[stride * 1] = Real(0.50979557910415916894194L) * (buffer[bufferStride * 1] + buffer[bufferStride * 0]);
    data[stride * 7] = Real(2.562915447741506178796086L) * (buffer[bufferStride * 1] - buffer[bufferStride * 0]);
}

template <MatrixOrScalarType T>
void dct_9(T* data, T* buffer, int stride, int bufferStride)
{
    buffer[bufferStride * 0] = data[stride * 6] - data[stride * 2];
    buffer[bufferStride * 1] = data[stride * 6] + data[stride * 2];
    buffer[bufferStride * 2] = data[stride * 7] - data[stride * 1];
    buffer[bufferStride * 3] = data[stride * 7] + data[stride * 1];

    buffer[bufferStride * 4] = data[stride * 0] - data[stride * 8];
    buffer[bufferStride * 5] = data[stride * 0] + data[stride * 8];
    buffer[bufferStride * 6] = data[stride * 5] - data[stride * 3];
    buffer[bufferStride * 7] = data[stride * 5] + data[stride * 3];

    data[stride * 2]         = buffer[bufferStride * 4] - buffer[bufferStride * 6];
    data[stride * 5]         = buffer[bufferStride * 6] + buffer[bufferStride * 4];
    buffer[bufferStride * 4] = buffer[bufferStride * 5] - buffer[bufferStride * 7];
    buffer[bufferStride * 6] = buffer[bufferStride * 7] + buffer[bufferStride * 5];

    data[stride * 3] =
        Real(0.866025403784438646763723170752936183471402627L) * (buffer[bufferStride * 0] + data[stride * 5]);
    buffer[bufferStride * 5] = data[stride * 4] + buffer[bufferStride * 3];
    buffer[bufferStride * 7] = buffer[bufferStride * 1] + buffer[bufferStride * 6];
    data[stride * 0]         = buffer[bufferStride * 5] + buffer[bufferStride * 7];
    data[stride * 6]         = Real(0.5) * FNMS(Real(2), buffer[bufferStride * 5], buffer[bufferStride * 7]);

    buffer[bufferStride * 5] =
        FNMS(Real(0.173648177666930348851716626769314796000375677L), data[stride * 2], buffer[bufferStride * 2]);
    buffer[bufferStride * 7] = FNMS(Real(0.5), data[stride * 5], buffer[bufferStride * 0]);
    buffer[bufferStride * 8] =
        FMA(Real(0.420276625461206169731530603237061658838781920L), buffer[bufferStride * 7], data[stride * 2]);
    data[stride * 5] =
        FNMS(Real(0.968908795874236621082202410917456709164223497L), buffer[bufferStride * 7], data[stride * 2]);
    data[stride * 1] =
        -(Real(0.866025403784438646763723170752936183471402627L) *
          (FNMS(Real(0.766044443118978035202392650555416673935832457L), data[stride * 5], buffer[bufferStride * 2])));
    data[stride * 5] =
        Real(0.866025403784438646763723170752936183471402627L) *
        (FMA(
            Real(0.939692620785908384054109277324731469936208134L), buffer[bufferStride * 8], buffer[bufferStride * 2]));
    data[stride * 7] =
        -(Real(0.984807753012208059366743024589523013670643252L) *
          (FNMS(
              Real(0.879385241571816768108218554649462939872416269L), buffer[bufferStride * 5],
              buffer[bufferStride * 7])));

    buffer[bufferStride * 0] = Real(0.5) * FNMS(Real(2), data[stride * 4], buffer[bufferStride * 3]);
    buffer[bufferStride * 7] = FNMS(Real(0.5), buffer[bufferStride * 6], buffer[bufferStride * 1]);
    buffer[bufferStride * 2] =
        FMA(Real(0.726681596905677465811651808188092531873167623L), buffer[bufferStride * 4], buffer[bufferStride * 7]);
    buffer[bufferStride * 5] =
        FNMS(Real(0.315207469095904627298647952427796244129086440L), buffer[bufferStride * 4], buffer[bufferStride * 7]);
    buffer[bufferStride * 6] =
        FNMS(Real(0.203604859554852403062088995281827210665664861L), buffer[bufferStride * 7], buffer[bufferStride * 4]);
    data[stride * 8] =
        FMS(Real(0.766044443118978035202392650555416673935832457L), buffer[bufferStride * 2], buffer[bufferStride * 0]);
    data[stride * 4] = -(
        FMA(Real(0.939692620785908384054109277324731469936208135L), buffer[bufferStride * 5], buffer[bufferStride * 0]));
    data[stride * 2] =
        FMA(Real(0.852868531952443209628250963940074071936020296L), buffer[bufferStride * 6], buffer[bufferStride * 0]);
}

template <MatrixOrScalarType T>
void dct_16(T* data, T* buffer, int stride, int bufferStride)
{
    buffer[bufferStride * 0] = data[stride * 0] + data[stride * 15];
    buffer[bufferStride * 1] = data[stride * 0] - data[stride * 15];
    buffer[bufferStride * 2] = data[stride * 8] - data[stride * 7];
    buffer[bufferStride * 3] = data[stride * 8] + data[stride * 7];
    data[stride * 8]         = buffer[bufferStride * 0] + buffer[bufferStride * 3];
    data[stride * 15]        = buffer[bufferStride * 0] - buffer[bufferStride * 3];

    buffer[bufferStride * 4] = data[stride * 6] - data[stride * 9];
    buffer[bufferStride * 5] = data[stride * 6] + data[stride * 9];
    buffer[bufferStride * 6] = data[stride * 1] - data[stride * 14];
    buffer[bufferStride * 7] = data[stride * 1] + data[stride * 14];

    buffer[bufferStride * 8] =
        FMA(Real(0.923879532511286756128183189396788286822416626L), buffer[bufferStride * 6],
            Real(0.382683432365089771728459984030398866761344562L) * buffer[bufferStride * 4]);
    buffer[bufferStride * 9] = FNMS(
        Real(0.382683432365089771728459984030398866761344562L), buffer[bufferStride * 6],
        Real(0.923879532511286756128183189396788286822416626L) * buffer[bufferStride * 4]);
    data[stride * 6] = buffer[bufferStride * 7] + buffer[bufferStride * 5];
    data[stride * 1] = buffer[bufferStride * 7] - buffer[bufferStride * 5];

    buffer[bufferStride * 4] = data[stride * 4] - data[stride * 11];
    buffer[bufferStride * 5] = data[stride * 4] + data[stride * 11];
    buffer[bufferStride * 6] = data[stride * 3] - data[stride * 12];
    buffer[bufferStride * 7] = data[stride * 3] + data[stride * 12];

    buffer[bufferStride * 10] =
        Real(0.707106781186547524400844362104849039284835938L) * (buffer[bufferStride * 4] + buffer[bufferStride * 6]);
    buffer[bufferStride * 11] =
        Real(0.707106781186547524400844362104849039284835938L) * (buffer[bufferStride * 4] - buffer[bufferStride * 6]);
    data[stride * 0] = buffer[bufferStride * 5] + buffer[bufferStride * 7];
    data[stride * 7] = buffer[bufferStride * 5] - buffer[bufferStride * 7];

    buffer[bufferStride * 12] = data[stride * 2] - data[stride * 13];
    buffer[bufferStride * 13] = data[stride * 2] + data[stride * 13];
    buffer[bufferStride * 0]  = data[stride * 10] - data[stride * 5];
    buffer[bufferStride * 3]  = data[stride * 10] + data[stride * 5];

    buffer[bufferStride * 14] = FNMS(
        Real(0.382683432365089771728459984030398866761344562L), buffer[bufferStride * 0],
        Real(0.923879532511286756128183189396788286822416626L) * buffer[bufferStride * 12]);
    buffer[bufferStride * 15] =
        FMA(Real(0.382683432365089771728459984030398866761344562L), buffer[bufferStride * 12],
            Real(0.923879532511286756128183189396788286822416626L) * buffer[bufferStride * 0]);
    data[stride * 14]         = buffer[bufferStride * 13] - buffer[bufferStride * 3];
    buffer[bufferStride * 3] += buffer[bufferStride * 13];

    buffer[bufferStride * 7]  = data[stride * 8] - data[stride * 0];
    buffer[bufferStride * 12] = buffer[bufferStride * 3] - data[stride * 6];
    data[stride * 4]          = FNMS(
        Real(0.382683432365089771728459984030398866761344563L), buffer[bufferStride * 12],
        Real(0.923879532511286756128183189396788286822416626L) * buffer[bufferStride * 7]);
    data[stride * 12] =
        FMA(Real(0.382683432365089771728459984030398866761344563L), buffer[bufferStride * 7],
            Real(0.923879532511286756128183189396788286822416626L) * buffer[bufferStride * 12]);
    buffer[bufferStride * 13] = data[stride * 8] + data[stride * 0];
    buffer[bufferStride * 0]  = buffer[bufferStride * 3] + data[stride * 6];
    data[stride * 8] =
        Real(0.707106781186547524400844362104849039284835938L) * (buffer[bufferStride * 13] - buffer[bufferStride * 0]);
    data[stride * 0] = buffer[bufferStride * 13] + buffer[bufferStride * 0];

    buffer[bufferStride * 0] =
        Real(0.707106781186547524400844362104849039284835938L) * (data[stride * 14] + data[stride * 1]);
    buffer[bufferStride * 6] = data[stride * 15] - buffer[bufferStride * 0];
    buffer[bufferStride * 7] = data[stride * 15] + buffer[bufferStride * 0];
    buffer[bufferStride * 3] =
        Real(0.707106781186547524400844362104849039284835938L) * (data[stride * 14] - data[stride * 1]);
    buffer[bufferStride * 12] = buffer[bufferStride * 3] - data[stride * 7];
    buffer[bufferStride * 13] = data[stride * 7] + buffer[bufferStride * 3];
    data[stride * 6]          = FNMS(
        Real(0.555570233019602224742830813948532874374937191L), buffer[bufferStride * 12],
        Real(0.831469612302545237078788377617905756738560812L) * buffer[bufferStride * 6]);
    data[stride * 14] =
        FMA(Real(0.980785280403230449126182236134239036973933731L), buffer[bufferStride * 13],
            Real(0.195090322016128267848284868477022240927691618L) * buffer[bufferStride * 7]);
    data[stride * 10] =
        FMA(Real(0.831469612302545237078788377617905756738560812L), buffer[bufferStride * 12],
            Real(0.555570233019602224742830813948532874374937191L) * buffer[bufferStride * 6]);
    data[stride * 2] = FNMS(
        Real(0.195090322016128267848284868477022240927691618L), buffer[bufferStride * 13],
        Real(0.980785280403230449126182236134239036973933731L) * buffer[bufferStride * 7]);

    buffer[bufferStride * 4]  = buffer[bufferStride * 1] - buffer[bufferStride * 10];
    buffer[bufferStride * 5]  = buffer[bufferStride * 15] - buffer[bufferStride * 9];
    buffer[bufferStride * 12] = buffer[bufferStride * 4] - buffer[bufferStride * 5];
    buffer[bufferStride * 13] = buffer[bufferStride * 4] + buffer[bufferStride * 5];
    buffer[bufferStride * 6]  = buffer[bufferStride * 14] - buffer[bufferStride * 8];
    buffer[bufferStride * 7]  = buffer[bufferStride * 11] - buffer[bufferStride * 2];
    buffer[bufferStride * 0]  = buffer[bufferStride * 6] - buffer[bufferStride * 7];
    buffer[bufferStride * 3]  = buffer[bufferStride * 7] + buffer[bufferStride * 6];

    data[stride * 5] = FNMS(
        Real(0.471396736825997648556387625905254377657460319L), buffer[bufferStride * 0],
        Real(0.881921264348355029712756863660388349508442620L) * buffer[bufferStride * 12]);
    data[stride * 13] =
        FMA(Real(0.956940335732208864935797886980269969482849206L), buffer[bufferStride * 3],
            Real(0.290284677254462367636192375817395274691476279L) * buffer[bufferStride * 13]);
    data[stride * 11] =
        FMA(Real(0.881921264348355029712756863660388349508442620L), buffer[bufferStride * 0],
            Real(0.471396736825997648556387625905254377657460319L) * buffer[bufferStride * 12]);
    data[stride * 3] = FNMS(
        Real(0.290284677254462367636192375817395274691476279L), buffer[bufferStride * 3],
        Real(0.956940335732208864935797886980269969482849206L) * buffer[bufferStride * 13]);

    buffer[bufferStride * 4]  = buffer[bufferStride * 1] + buffer[bufferStride * 10];
    buffer[bufferStride * 5]  = buffer[bufferStride * 14] + buffer[bufferStride * 8];
    buffer[bufferStride * 12] = buffer[bufferStride * 4] + buffer[bufferStride * 5];
    buffer[bufferStride * 3]  = buffer[bufferStride * 4] - buffer[bufferStride * 5];
    buffer[bufferStride * 7]  = buffer[bufferStride * 2] + buffer[bufferStride * 11];
    buffer[bufferStride * 6]  = buffer[bufferStride * 15] + buffer[bufferStride * 9];
    buffer[bufferStride * 13] = buffer[bufferStride * 7] + buffer[bufferStride * 6];
    buffer[bufferStride * 0]  = buffer[bufferStride * 6] - buffer[bufferStride * 7];

    data[stride * 1] = FNMS(
        Real(0.098017140329560601994195563888641845861136673L), buffer[bufferStride * 13],
        Real(0.995184726672196886244836953109479921575474869L) * buffer[bufferStride * 12]);
    data[stride * 9] =
        FMA(Real(0.634393284163645498215171613225493370675687095L), buffer[bufferStride * 3],
            Real(0.773010453362736960810906609758469800971041293L) * buffer[bufferStride * 0]);
    data[stride * 15] =
        FMA(Real(0.098017140329560601994195563888641845861136673L), buffer[bufferStride * 12],
            Real(0.995184726672196886244836953109479921575474869L) * buffer[bufferStride * 13]);
    data[stride * 7] = FNMS(
        Real(0.634393284163645498215171613225493370675687095L), buffer[bufferStride * 0],
        Real(0.773010453362736960810906609758469800971041293L) * buffer[bufferStride * 3]);
}

#undef FMA
#undef FMS
#undef FNMA
#undef FNMS

// buffer must have size bufferStride*N
template <MatrixOrScalarType T>
void dct(T* x, int N, T* buffer, int stride, int bufferStride)
{
    switch (N)
    {
    case 0:
        return;
    case 1:
        return;
    case 3:
        dct_3(x, buffer, stride, bufferStride);
        return;
    case 6:
        dct_6(x, buffer, stride, bufferStride);
        return;
    case 8:
        dct_8(x, buffer, stride, bufferStride);
        return;
    case 9:
        dct_9(x, buffer, stride, bufferStride);
        return;
    case 16:
        dct_16(x, buffer, stride, bufferStride);
        return;
    }

    // Check that N is power of 2
    if (N % 2 == 1)
    {
        throw std::runtime_error("dct only implemented for N=2^m, N=3*2^m, N=9*2^m.");
    }

    int halfLen = N / 2;

    for (int i = 0; i < halfLen; ++i)
    {
        buffer[bufferStride * i] = x[stride * i] + x[stride * (N - 1 - i)];
        buffer[bufferStride * (i + halfLen)] =
            (x[stride * i] - x[stride * (N - 1 - i)]) / (std::cos(Real(i + 0.5) * std::numbers::pi_v<Real> / N) * 2);
    }

    dct(buffer, halfLen, x, bufferStride, stride);
    dct(&buffer[bufferStride * halfLen], halfLen, x, bufferStride, stride);

    for (int i = 0; i < halfLen - 1; ++i)
    {
        x[stride * (i * 2 + 0)] = buffer[bufferStride * i];
        x[stride * (i * 2 + 1)] = buffer[bufferStride * (i + halfLen)] + buffer[bufferStride * (i + halfLen + 1)];
    }

    x[stride * (N - 2)] = buffer[bufferStride * (halfLen - 1)];
    x[stride * (N - 1)] = buffer[bufferStride * (N - 1)];
}

///
/// \brief      Computes the discrete cosine transform (type II) of _x_.
///
/// \copybrief dct(T* x, int N, T* buffer)
/// This is defined as
///
/// \f[
///     X_k = \sum_{n=0}^{N-1} x_n \cos \left[ \frac{\pi}{N} \left(n + \frac{1}{2} \right) k \right].
/// \f]
///
/// This implementation does not work in-place, but instead needs a buffer of the same length as the input. For small
/// \f$N\f$ there are specifically optimized codelets. For large \f$N\f$ a recursive approach is used inspired by
/// https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms. The function is currently only implemented for
/// \f$N=k\cdot 2^m\f$ with \f$k\in\{1,3,9\}\f$ and \f$m\geq 0\f$.
///
/// \param      x       The data that is transformed.
/// \param      N       The number of elements to which the DCT is applied.
/// \param      buffer  Buffer of at least length \a bufferStride*N.
///
template <MatrixOrScalarType T>
void dct(T* x, int N, T* buffer)
{
    dct(x, N, buffer, 1, 1);
}

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
extern template void dct<Real>(Real*, int, Real*);
extern template void dct<Complex>(Complex*, int, Complex*);
extern template void dct<RealVector>(RealVector*, int, RealVector*);
extern template void dct<Vector>(Vector*, int, Vector*);
extern template void dct<RealMatrix>(RealMatrix*, int, RealMatrix*);
extern template void dct<Matrix>(Matrix*, int, Matrix*);
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

///
/// \brief      Convenience version that automatically allocates the required buffer.
///
template <MatrixOrScalarType T>
void dct(T* x, int N)
{
    std::vector<T> buffer(N);
    dct(x, N, buffer.data(), 1, 1);
}

// clang-format off
#ifndef SCICORE_DONT_PRECOMPILE_TEMPLATES
extern template void dct<Real>(Real*, int);
extern template void dct<Complex>(Complex*, int);
extern template void dct<RealVector>(RealVector*, int);
extern template void dct<Vector>(Vector*, int);
extern template void dct<RealMatrix>(RealMatrix*, int);
extern template void dct<Matrix>(Matrix*, int);
#endif // SCICORE_DONT_PRECOMPILE_TEMPLATES
// clang-format on

/// \} // end of DCT

} // namespace SciCore

#endif // SCICORE_DCT_H
