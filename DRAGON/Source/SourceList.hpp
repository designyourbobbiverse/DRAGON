//
//  SourceList.hpp
//  DRAGON/Source
//
//  Created by Bobbie Markwick on 05/09/2026.
//

#ifdef SOURCE_LIST_HPP //Guard against accidental include

#include <vector>

#include <memory> //std::unique_ptr
#include <utility> //std::forward
#include <concepts> //std::derived_from
#include <type_traits> //std::decay_t

namespace DRAGON::Source {
class SourceList;
template<class S> concept SourceType = std::derived_from<std::decay_t<S>, SourceTerm> && !std::derived_from<std::decay_t<S>, SourceList>;


class SourceList : public SourceTerm {
private:
    std::vector<std::unique_ptr<SourceTerm>>  sources;
public:
    ConservativeState source_density(const PrimitiveState &w, double t) override {
        ConservativeState S{};
        for(auto& s: sources) S += s->source_density(w,t);
        return S;
    }
    
    //MARK: Construction
    template<SourceType... Ss> SourceList(Ss&&... ss) {
        (append(std::forward<Ss>(ss)), ...);
    }
    
    SourceList(const SourceList&) = delete; //No copying
    SourceList(SourceList&&) noexcept = default;
    SourceList& operator=(SourceList&&) noexcept = default;
    template<SourceType S> SourceList& operator=(S&& s){
        clear();
        append(std::forward<S>(s));
        return *this;
    }
    
    
    //MARK: Append
    //Add an element to the end of the list
    template<SourceType S> void append(S&& s){
        static_assert(std::derived_from<std::decay_t<S>, SourceTerm>);
        sources.push_back(std::make_unique<std::decay_t<S>>(std::forward<S>(s)));
    }
    //Move elements from another list to the end of this list
    inline void append(SourceList&& ss){
        for (auto& s : ss.sources)  sources.push_back(std::move(s));
        ss.clear();
    }
    
    //MARK: Prepend
    //Add an element to the beginning of the list
    template<SourceType S> void prepend(S&& s){
        static_assert(std::derived_from<std::decay_t<S>, SourceTerm>);
        sources.insert(sources.begin(), std::make_unique<std::decay_t<S>>(std::forward<S>(s)));
    }
    //Move elements from another list to the beginning of this list
    inline void prepend(SourceList&& ss){
        ss.append(std::move(*this));
        sources = std::move(ss.sources);
    }
    //MARK: Clear
    //Remove all elements from the list
    void clear(){ sources.clear(); }
    
};
}


//MARK: + Operators
template<DRAGON::Source::SourceType A, DRAGON::Source::SourceType B>  DRAGON::Source::SourceList operator+(A&& a, B&& b){
    DRAGON::Source::SourceList result;
    result.append(std::forward<A>(a));
    result.append(std::forward<B>(b));
    return result;
}
template<DRAGON::Source::SourceType B> DRAGON::Source::SourceList& operator+=(DRAGON::Source::SourceList& lhs, B&& rhs){
    append(std::forward<B>(rhs));
    return lhs;
}
template<DRAGON::Source::SourceType B> DRAGON::Source::SourceList operator+(DRAGON::Source::SourceList lhs, B&& rhs){ lhs += rhs; return lhs; }
template<DRAGON::Source::SourceType A> DRAGON::Source::SourceList operator+(A&& lhs, DRAGON::Source::SourceList rhs){
    rhs.prepend(std::forward<A>(lhs));
    return rhs;
}
inline DRAGON::Source::SourceList& operator+=(DRAGON::Source::SourceList& lhs, DRAGON::Source::SourceList rhs){
    lhs.append(std::move(rhs));
    return lhs;
}
inline DRAGON::Source::SourceList operator+(DRAGON::Source::SourceList lhs, DRAGON::Source::SourceList rhs) {
    lhs += std::move(rhs);
    return lhs;
}

#else
#warning SourceList.hpp appears to have been included accidentally. Typically, you want to include Source.hpp instead.
#endif
