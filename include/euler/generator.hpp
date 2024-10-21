#pragma once

#include <coroutine>
#include <ranges>
#include <utility>

inline namespace euler
{
template <typename T> class generator
{
    struct promise
    {
        promise() = default;

        generator get_return_object() { return generator(std::coroutine_handle<promise>::from_promise(*this)); }

        [[nodiscard]] std::suspend_always initial_suspend() const { return {}; }
        [[nodiscard]] std::suspend_always final_suspend() const noexcept { return {}; }

        void return_void() const noexcept {}

        void unhandled_exception() noexcept { exception_ = std::current_exception(); }

        void rethrow_if_exception()
        {
            if (exception_)
                std::rethrow_exception(exception_);
        }

        std::suspend_always yield_value(auto &&v) noexcept
        {
            value_ = v;
            return {};
        }

        std::exception_ptr exception_;
        T value_;
    };

  public:
    using promise_type = promise;
    class sentinel
    {
    };

    class iterator
    {
        using handle_type = std::coroutine_handle<promise_type>;

      public:
        using difference_type = std::ptrdiff_t;

        iterator() = default;
        ~iterator()
        {
            if (handle_)
                handle_.destroy();
        }

        // Non-copyable because coroutine handles point to a unique resource
        iterator(iterator const &) = delete;
        iterator(iterator &&rhs) noexcept : handle_(std::exchange(rhs.handle_, nullptr)) {}
        iterator &operator=(iterator const &) = delete;
        iterator &operator=(iterator &&rhs) noexcept
        {
            handle_ = std::exchange(rhs.handle_, nullptr);
            return *this;
        }

        friend bool operator==(iterator const &it, sentinel /*unused*/) noexcept
        {
            return (!it.handle_ || it.handle_.done());
        }

        iterator &operator++()
        {
            handle_.resume();
            if (handle_.done())
                handle_.promise().rethrow_if_exception();
            return *this;
        }

        void operator++(int) { (void)this->operator++(); }

        T &operator*() const { return handle_.promise().value_; }

      private:
        friend class generator;
        explicit iterator(handle_type handle) : handle_(handle) {}

        handle_type handle_;
    };

    using handle_type = std::coroutine_handle<promise_type>;

    generator() noexcept = default;
    ~generator()
    {
        if (handle_)
            handle_.destroy();
    }

    generator(generator const &) = delete;
    generator(generator &&rhs) noexcept : handle_(std::exchange(rhs.handle_, nullptr)) {}
    generator &operator=(generator const &) = delete;
    generator &operator=(generator &&rhs) noexcept
    {
        swap(rhs);
        return *this;
    }

    iterator begin()
    {
        handle_.resume();
        if (handle_.done())
            handle_.promise().rethrow_if_exception();
        return {std::exchange(handle_, nullptr)};
    }

    sentinel end() const noexcept { return {}; }

    void swap(generator &other) noexcept { std::swap(handle_, other.handle_); }

  private:
    friend class iterator;
    explicit generator(handle_type handle) noexcept : handle_(handle) {}

    handle_type handle_ = nullptr;
};
} // namespace euler

template <typename T> constexpr bool std::ranges::enable_view<generator<T>> = true;
