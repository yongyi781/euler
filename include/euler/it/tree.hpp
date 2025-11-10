#pragma once

#include "base.hpp"

namespace euler
{
namespace it
{
/// Predicate that always returns true.
constexpr auto pred_true = [](auto &&) { return true; };

/// Enumerates a virtual tree. This version enumerates iteratively, rather than recursively.
/// Usage: `it::tree(root, fun(x, f, <call f on each child of x>), fun(x, <whether x has children>))`.
template <typename T, std::invocable<T, result_t(T)> Fun, std::predicate<T> Pred = decltype(pred_true)>
class tree : public it_base
{
    T _root;
    Fun _childrenFun;
    Pred _hasChildren;

  public:
    using value_type = T;

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    /// @param hasChildren A function that returns whether a node has children. Optional, but really speeds up the
    ///                    runtime if provided.
    constexpr tree(T root, Fun childrenFun, Pred hasChildren)
        : _root(std::move(root)), _childrenFun(std::move(childrenFun)), _hasChildren(std::move(hasChildren))
    {
    }

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    constexpr tree(T root, Fun childrenFun) : tree(std::move(root), std::move(childrenFun), pred_true) {}

    template <typename Callback> constexpr result_t operator()(Callback f) const
    {
        if (!it::callbackResult(f, _root))
            return result_break;
        std::vector<T> s{std::move(_root)};
        while (!s.empty())
        {
            T current = std::move(s.back());
            s.pop_back();
            if (!it::callbackResult(_childrenFun, std::move(current), [&](T child) -> result_t {
                    if (!it::callbackResult(f, child))
                        return result_break;
                    if (std::invoke(_hasChildren, child))
                        s.emplace_back(std::move(child));
                    return result_continue;
                }))
                return result_break;
        }
        return result_continue;
    }
};

/// Enumerates a virtual tree. This version enumerates recursively and in preorder traversal order.
/// Usage: `it::tree_preorder(root, fun(x, f, <call f on each child of x>), fun(x, <whether x has children>))`.
template <typename T, std::invocable<T, result_t(T)> Fun, std::predicate<T> Pred = decltype(pred_true)>
class tree_preorder : public it_base
{
    T _root;
    Fun _childrenFun;
    Pred _hasChildren;

    template <std::invocable<T> Callback> constexpr result_t dfs(T node, Callback f) const
    {
        return callbackResult(_childrenFun, std::move(node), [&](T child) -> result_t {
            if (!callbackResult(f, child))
                return result_break;
            if (!std::invoke(_hasChildren, child))
                return result_continue;
            return dfs(std::move(child), f);
        });
    }

  public:
    using value_type = T;

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    /// @param hasChildren A function that returns whether a node has children. Optional, but really speeds up the
    ///                    runtime if provided.
    constexpr tree_preorder(T root, Fun childrenFun, Pred hasChildren)
        : _root(std::move(root)), _childrenFun(std::move(childrenFun)), _hasChildren(std::move(hasChildren))
    {
    }

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    constexpr tree_preorder(T root, Fun childrenFun) : tree_preorder(std::move(root), std::move(childrenFun), pred_true)
    {
    }

    template <typename Callback> constexpr result_t operator()(Callback f) const
    {
        if (!callbackResult(f, _root))
            return result_break;
        return dfs(_root, std::move(f));
    }
};
} // namespace it
} // namespace euler
