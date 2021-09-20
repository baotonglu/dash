#include "dash.h"

extern "C" tree_api* create_tree(const tree_options_t& opt)
{
    return new dash_wrapper();
}