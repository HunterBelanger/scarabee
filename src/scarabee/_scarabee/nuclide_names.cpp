#include <utils/nuclide_names.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <cctype>

namespace scarabee {

std::string nuclide_name_to_internal_name(const std::string& name) {
  // Must remove any _ from nuclide name (for TSLs like H1_H2O)
  // Also discards isomer information (Ag110m1 -> Ag110)

  std::string simp_name;
  simp_name.reserve(name.size());

  bool doing_nums = false;
  for (std::size_t i = 0; i < 5; i++) {
    const bool i_is_digit = std::isdigit(name[i]);

    if (i_is_digit) doing_nums = true;
    if (doing_nums && i_is_digit == false) break;

    simp_name += name[i];
  }

  return simp_name;
}

std::string nuclide_name_to_simple_name(const std::string& name) {
  // Must remove any _ from nuclide name (for TSLs like H1_H2O)
  std::string simp_name = name;
  auto loc_undr_scr = simp_name.find('_');
  if (loc_undr_scr != std::string::npos) {
    simp_name.resize(loc_undr_scr);
  }
  return simp_name;
}

std::string nuclide_name_to_element_symbol(const std::string& name) {
  std::string elem_name;

  for (std::size_t i = 0; i < 2; i++) {
    if (std::isalpha(name[i])) {
      elem_name += name[i];
    } else {
      break;
    }
  }

  return elem_name;
}

std::uint32_t nuclide_name_to_za(const std::string& name) {
  const char first_letter = name[0];
  const char second_letter = name[1];

  const bool second_is_letter = std::isalpha(second_letter);

  // First, find the element
  std::uint32_t Z = 1;
  bool found_element = false;
  for (Z = 1; Z < ELEMENTS.size(); Z++) {
    const auto& Elem = ELEMENTS[Z];

    if (first_letter == Elem.symbol[0]) {
      // First letters match ! What about second ?
      if ((second_is_letter && second_letter == Elem.symbol[1]) ||
          (second_is_letter == false)) {
        found_element = true;
        break;
      }
    }
  }

  if (found_element == false) {
    const auto mssg = "Could not find an element for \"" + name + "\".";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Now we need to get the A
  std::uint32_t A = 0;
  std::uint32_t m = 0;
  bool store_A = false;
  bool store_m = false;
  for (const auto c : name) {
    const bool c_is_digit = std::isdigit(c);

    if (store_m && c_is_digit) {
      m = static_cast<std::uint32_t>(c - '0');

      // Should only be one digit ! We can now break
      break;
    } else if (store_m) {
      // Next character wasn't a digit. This is a poorly formed name.
      const auto mssg = "The name \"" + name + "\" is not a valid nuclide.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }

    if (store_A == false and c_is_digit) {
      store_A = true;
    }

    if (store_A && c_is_digit) {
      A *= 10;  // Move things over one decimal place. Doesn't matter for first
                // number, as A is 0 !
      A += static_cast<std::uint32_t>(c - '0');
    } else if (store_A) {
      // Apparently we have gone through all digits for A
      store_A = false;

      if (c == 'm') {
        // If c is m, then we have an isomer state to obtain !
        store_m = true;
      } else {
        // No isomer state, we can exit
        break;
      }
    }
  }

  return (Z * 1000 + A) * 10 + m;
}

}  // namespace scarabee