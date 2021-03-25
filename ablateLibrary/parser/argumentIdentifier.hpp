#ifndef ABLATELIBRARY_ARGUMENTIDENTIFIER_HPP
#define ABLATELIBRARY_ARGUMENTIDENTIFIER_HPP
#include <optional>
#include <string>

#define TMP_COMMA ,

#define ARG(interfaceTypeFullName, inputName, description) \
    ablate::parser::ArgumentIdentifier<interfaceTypeFullName> { inputName, description, false }

#define OPT(interfaceTypeFullName, inputName, description) \
    ablate::parser::ArgumentIdentifier<interfaceTypeFullName> { inputName, description, true }

namespace ablate::parser {
template <typename Interface>
struct ArgumentIdentifier {
    const std::string inputName;
    const std::string description;
    const bool optional;
    bool operator==(const ArgumentIdentifier<Interface>& other) const { return inputName == other.inputName && description == other.description && optional == other.optional; }
};
}  // namespace ablate::parser

#endif  // ABLATELIBRARY_ARGUMENTIDENTIFIER_HPP
