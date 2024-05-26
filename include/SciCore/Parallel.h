//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#ifndef SCICORE_PARALLEL_H
#define SCICORE_PARALLEL_H

#include <concepts>
#include <limits>

#include <taskflow/taskflow.hpp>

#include "Definitions.h"
#include "SciCore_export.h"

namespace SciCore
{

template <typename FunctionT>
    requires std::invocable<FunctionT, int>
SCICORE_EXPORT void parallelFor(
    FunctionT&& f,
    int begin,
    int end,
    tf::Executor& executor,
    int nChunks = std::numeric_limits<int>::max())
{
    if (nChunks <= 1 || executor.num_workers() <= 1)
    {
        for (int j = begin; j < end; ++j)
        {
            f(j);
        }
        return;
    }

    if (nChunks > end - begin)
    {
        nChunks = end - begin;
    }

    tf::Taskflow taskflow;

    int chunkSize = (end - begin) / nChunks;
    for (int i = 0; i < nChunks - 1; ++i)
    {
        taskflow.emplace(
            [&f, begin, chunkSize, i]()
            {
                int currentBegin = begin + i * chunkSize;
                int currentEnd   = begin + (i + 1) * chunkSize;
                for (int j = currentBegin; j < currentEnd; ++j)
                {
                    f(j);
                }
            });
    }

    taskflow.emplace(
        [&f, begin, end, nChunks, chunkSize]()
        {
            int currentBegin = begin + (nChunks - 1) * chunkSize;
            for (int j = currentBegin; j < end; ++j)
            {
                f(j);
            }
        });

    executor.run(taskflow).get();
}

template <typename FunctionT>
    requires std::invocable<FunctionT, int>
SCICORE_EXPORT auto parallelSum(
    FunctionT&& f,
    int begin,
    int end,
    tf::Executor& executor,
    int nChunks = std::numeric_limits<int>::max())
{
#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    // 64 bytes on x86-64
    constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

    using ReturnType = std::decay_t<typename std::invoke_result<FunctionT, int>::type>;

    if (nChunks <= 1 || executor.num_workers() <= 1)
    {
        ReturnType returnValue = f(begin);
        for (int j = begin + 1; j < end; ++j)
        {
            returnValue += f(j);
        }
        return returnValue;
    }

    if (nChunks > end - begin)
    {
        nChunks = end - begin;
    }

    tf::Taskflow taskflow;

    struct ResultStruct
    {
        alignas(hardware_destructive_interference_size)
            ReturnType value; // FIXME necessary? supposed to prevent false sharing
    };

    std::vector<ResultStruct> results(nChunks);

    int chunkSize = (end - begin) / nChunks;
    for (int i = 0; i < nChunks - 1; ++i)
    {
        taskflow.emplace(
            [&f, &results, begin, chunkSize, i]()
            {
                int currentBegin         = begin + i * chunkSize;
                int currentEnd           = begin + (i + 1) * chunkSize;
                ReturnType partialResult = f(currentBegin);
                for (int j = currentBegin + 1; j < currentEnd; ++j)
                {
                    partialResult += f(j);
                }
                results[i].value = std::move(partialResult);
            });
    }

    taskflow.emplace(
        [&f, &results, begin, end, nChunks, chunkSize]()
        {
            int currentBegin         = begin + (nChunks - 1) * chunkSize;
            ReturnType partialResult = f(currentBegin);
            for (int j = currentBegin + 1; j < end; ++j)
            {
                partialResult += f(j);
            }
            results[nChunks - 1].value = std::move(partialResult);
        });

    executor.run(taskflow).get();
    ReturnType returnValue = std::move(results[0].value);
    for (size_t i = 1; i < results.size(); ++i)
    {
        returnValue += results[i].value;
    }
    return returnValue;
}

} // namespace SciCore

#endif // SCICORE_PARALLEL_H
