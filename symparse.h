#pragma once
// symparse.h  –  version_2026
//
// Polynomial expression pre-processor for the right-hand-side parser.
//
// The parser handles most polynomial forms natively, including
// powers  x^2,  x*y,  (x+y)^2,  (x-y)^3,  etc.  This file therefore
// provides only lightweight string sanitisation (whitespace removal, basic
// implicit-multiplication insertion).
//
// For heavy symbolic expansion (e.g. (x+y)^20) compile with:
//     -DHAVE_SYMENGINE   and link against SymEngine.
//
// Without HAVE_SYMENGINE the function performs basic clean-up, which is
// sufficient for all built-in examples.

#include <string>
#include <cctype>
#include <stdexcept>

#ifdef HAVE_SYMENGINE
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/symbol.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/parser.h>
using SymEngine::expand;
#endif

// ───────────────────────────────────────────────────────────────────────────
//  Lightweight version (no SymEngine)
// ───────────────────────────────────────────────────────────────────────────

// Remove whitespace and insert '*' between a digit/letter and a following
// letter (implicit multiplication like "2x" → "2*x").
inline std::string basic_sanitize(const std::string& input) {
    std::string s;
    s.reserve(input.size() * 2);
    for (size_t i = 0; i < input.size(); ++i) {
        char c = input[i];
        if (c == ' ' || c == '\t' || c == '\n') continue;
        // Insert '*' when a digit or ')' is immediately followed by a letter or '('
        if (!s.empty()) {
            char prev = s.back();
            bool digitThenAlpha = std::isdigit(static_cast<unsigned char>(prev)) &&
                                  std::isalpha(static_cast<unsigned char>(c));
            bool parenThenAlnum = (prev == ')') &&
                                  (std::isalpha(static_cast<unsigned char>(c)) ||
                                   std::isdigit(static_cast<unsigned char>(c)) ||
                                   c == '(');
            if (digitThenAlpha || parenThenAlnum)
                s += '*';
        }
        s += c;
    }
    return s;
}

// ───────────────────────────────────────────────────────────────────────────
//  Public API
// ───────────────────────────────────────────────────────────────────────────

// Pre-process a polynomial expression string for the right-hand-side parser.
// - Without HAVE_SYMENGINE: whitespace cleanup + implicit multiplication.
// - With    HAVE_SYMENGINE: full algebraic expansion via SymEngine.
inline std::string expand_expression(const std::string& input) {
#ifdef HAVE_SYMENGINE
    // Replace ^ with ** for SymEngine
    std::string adjusted = input;
    {
        size_t pos = 0;
        while ((pos = adjusted.find('^', pos)) != std::string::npos) {
            adjusted.replace(pos, 1, "**");
            pos += 2;
        }
        // Remove spaces
        std::string tmp;
        for (char c : adjusted) if (c != ' ') tmp += c;
        adjusted = tmp;
    }
    try {
        auto expr     = SymEngine::parse(adjusted);
        auto expanded = expand(expr);
        std::string result = expanded->__str__();
        // Replace ** back to ^
        {
            size_t pos = 0;
            while ((pos = result.find("**", pos)) != std::string::npos) {
                result.replace(pos, 2, "^");
                pos += 1;
            }
        }
        // Insert * between digit and letter
        return basic_sanitize(result);
    } catch (const std::exception& e) {
        // Fall back to basic sanitise
        return basic_sanitize(input);
    }
#else
    return basic_sanitize(input);
#endif
}
