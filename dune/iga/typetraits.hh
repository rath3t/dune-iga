//
// Created by lex on 14/12/2021.
//

#pragma once

template < template <typename...> class Template, typename T >
struct is_instantiation_of : std::false_type {};

template < template <typename...> class Template, typename... Args >
struct is_instantiation_of< Template, Template<Args...> > : std::true_type {};