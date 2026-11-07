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

#include <algorithm>
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
//  Constants in scientific notation
// ───────────────────────────────────────────────────────────────────────────

// The right-hand-side parser reads plain decimal constants only, so every
// constant written in scientific notation, m e k, is rewritten as
// (m*10^k) or (m/10^|k|); the power of ten is then evaluated in interval
// arithmetic.  Names of variables and functions are copied unchanged, so
// "x1e2" stays a name and "2e" or "2ex" are not read as constants.
inline std::string rewrite_scientific(const std::string& in) {
    auto isDigit = [](char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; };
    auto isName  = [](char c) { return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_'; };

    std::string out;
    out.reserve(in.size() + 16);
    const size_t N = in.size();
    size_t i = 0;
    while (i < N) {
        const char c = in[i];
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {      // a name
            size_t j = i;
            while (j < N && isName(in[j])) ++j;
            out.append(in, i, j - i);
            i = j;
            continue;
        }
        if (isDigit(c) || (c == '.' && i + 1 < N && isDigit(in[i + 1]))) {  // a constant
            size_t j = i;
            while (j < N && (isDigit(in[j]) || in[j] == '.')) ++j;
            const std::string mantissa = in.substr(i, j - i);
            if (j < N && (in[j] == 'e' || in[j] == 'E')) {
                size_t k = j + 1;
                char sign = '+';
                if (k < N && (in[k] == '+' || in[k] == '-')) sign = in[k++];
                if (k < N && isDigit(in[k])) {
                    size_t e0 = k;
                    while (k < N && isDigit(in[k])) ++k;
                    std::string expo = in.substr(e0, k - e0);
                    expo.erase(0, std::min(expo.find_first_not_of('0'), expo.size() - 1));
                    out += "(" + mantissa + (sign == '-' ? "/" : "*") + "10^" + expo + ")";
                    i = k;
                    continue;
                }
            }
            out += mantissa;
            i = j;
            continue;
        }
        out += c;
        ++i;
    }
    return out;
}

// ───────────────────────────────────────────────────────────────────────────
//  Lightweight version (no SymEngine)
// ───────────────────────────────────────────────────────────────────────────

// Remove whitespace and insert '*' for implicit multiplication: after a
// constant that is followed by a letter or '(' ("2x" → "2*x",
// "2(x+1)" → "2*(x+1)"), and after ')' that is followed by a letter, a digit
// or '(' ("(x+1)(y-1)" → "(x+1)*(y-1)").  Digits inside a name ("x1", "k2a",
// "atan2") are part of the name.
inline std::string basic_sanitize(const std::string& input) {
    // Does s end in a constant, i.e. in digits that do not belong to a name?
    auto endsInConstant = [](const std::string& s) {
        size_t k = s.size();
        while (k > 0 && (std::isdigit(static_cast<unsigned char>(s[k - 1])) || s[k - 1] == '.')) --k;
        if (k == s.size()) return false;
        return k == 0 || !(std::isalpha(static_cast<unsigned char>(s[k - 1])) || s[k - 1] == '_');
    };
    std::string s;
    s.reserve(input.size() * 2);
    for (size_t i = 0; i < input.size(); ++i) {
        char c = input[i];
        if (c == ' ' || c == '\t' || c == '\n') continue;
        // Insert '*' when a constant or ')' is immediately followed by a letter or '('
        if (!s.empty()) {
            char prev = s.back();
            bool digitThenAlpha = std::isdigit(static_cast<unsigned char>(prev)) &&
                                  (std::isalpha(static_cast<unsigned char>(c)) || c == '(') &&
                                  endsInConstant(s);
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
// - Both:                   scientific notation rewritten, see above.
// - Without HAVE_SYMENGINE: whitespace cleanup + implicit multiplication.
// - With    HAVE_SYMENGINE: full algebraic expansion via SymEngine.
inline std::string expand_expression(const std::string& raw) {
    const std::string input = rewrite_scientific(raw);
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
