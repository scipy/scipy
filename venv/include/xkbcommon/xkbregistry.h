/*
 * Copyright © 2020 Red Hat, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */


#ifndef _XKBREGISTRY_H_
#define _XKBREGISTRY_H_

#include <stdarg.h>
#include <stdbool.h>

/**
 * @file
 * @brief Query for available RMLVO
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup registry Query for available RMLVO
 *
 * The libxkbregistry API to query for available rules, models, layouts,
 * variants and options (RMLVO). libxkbregistry is a separate library to
 * libxkbcommon.
 *
 * This library is the replacement for clients currently parsing evdev.xml
 * directly. The library is intended to provide easy access to the set of
 * **possible** MLVO configurations for a given ruleset. It is not a library to
 * apply these configurations, merely to enumerate them. The intended users of
 * this library are the configuration UIs that allow a user to select their
 * keyboard layout of choice.
 *
 * @{
 */

/**
 * @struct rxkb_context
 *
 * Opaque top level library context object.
 *
 * The context contains general library state, like include paths and parsed
 * data. Objects are created in a specific context, and multiple contexts
 * may coexist simultaneously. Objects from different contexts are
 * completely separated and do not share any memory or state.
 */
struct rxkb_context;

/**
 * @struct rxkb_model
 *
 * Opaque struct representing an XKB model.
 */
struct rxkb_model;

/**
 * @struct rxkb_layout
 *
 * Opaque struct representing an XKB layout, including an optional variant.
 * Where the variant is NULL, the layout is the base layout.
 *
 * For example, "us" is the base layout, "us(intl)" is the "intl" variant of the
 * layout "us".
 */
struct rxkb_layout;

/**
 * @struct rxkb_option_group
 *
 * Opaque struct representing an option group. Option groups divide the
 * individual options into logical groups. Their main purpose is to indicate
 * whether some options are mutually exclusive or not.
 */
struct rxkb_option_group;

/**
 * @struct rxkb_option
 *
 * Opaque struct representing an XKB option. Options are grouped inside an @ref
 * rxkb_option_group.
 */
struct rxkb_option;

/**
 *
 * @struct rxkb_iso639_code
 *
 * Opaque struct representing an ISO 639-3 code (e.g. "eng", "fra"). There
 * is no guarantee that two identical ISO codes share the same struct. You
 * must not rely on the pointer value of this struct.
 *
 * See https://iso639-3.sil.org/code_tables/639/data for a list of codes.
 */
struct rxkb_iso639_code;

/**
 *
 * @struct rxkb_iso3166_code
 *
 * Opaque struct representing an ISO 3166 Alpha 2 code (e.g. "US", "FR").
 * There is no guarantee that two identical ISO codes share the same struct.
 * You must not rely on the pointer value of this struct.
 *
 * See https://en.wikipedia.org/wiki/List_of_ISO_3166_country_codes for a list
 * of codes.
 */
struct rxkb_iso3166_code;

/**
 * Describes the popularity of an item. Historically, some highly specialized or
 * experimental definitions are excluded from the default list and shipped in
 * separate files. If these extra definitions are loaded (see @ref
 * RXKB_CONTEXT_LOAD_EXOTIC_RULES), the popularity of the item is set
 * accordingly.
 *
 * If the exotic items are not loaded, all items will have the standard
 * popularity.
 */
enum rxkb_popularity {
    RXKB_POPULARITY_STANDARD = 1,
    RXKB_POPULARITY_EXOTIC,
};

/**
 * Flags for context creation.
 */
enum rxkb_context_flags {
    RXKB_CONTEXT_NO_FLAGS = 0,
    /**
     * Skip the default include paths. This requires the caller to call
     * rxkb_context_include_path_append() or
     * rxkb_context_include_path_append_default().
     */
    RXKB_CONTEXT_NO_DEFAULT_INCLUDES = (1 << 0),
    /**
     * Load the extra items that are considered too exotic for the default list.
     *
     * For historical reasons, xkeyboard-config ships those exotic rules in a
     * separate file (e.g. `evdev.extras.xml`). Where the exotic rules are
     * requested, libxkbregistry will look for and load `$ruleset.extras.xml`
     * in the include paths, see rxkb_context_include_path_append() for details
     * on the lookup behavior.
     */
    RXKB_CONTEXT_LOAD_EXOTIC_RULES = (1 << 1),
    /**
     * Disable the use of secure_getenv for this context, so that privileged
     * processes can use environment variables. Client uses at their own risk.
     *
     * @since 1.5.0
     */
    RXKB_CONTEXT_NO_SECURE_GETENV = (1 << 2)
};

/**
 * Create a new xkb registry context.
 *
 * The context has an initial refcount of 1. Use rxkb_context_unref() to release
 * memory associated with this context.
 *
 * Creating a context does not parse the files yet, use
 * rxkb_context_parse().
 *
 * @param flags Flags affecting context behavior
 * @return A new xkb registry context or NULL on failure
 */
struct rxkb_context *
rxkb_context_new(enum rxkb_context_flags flags);

/** Specifies a logging level. */
enum rxkb_log_level {
    RXKB_LOG_LEVEL_CRITICAL = 10, /**< Log critical internal errors only. */
    RXKB_LOG_LEVEL_ERROR = 20,    /**< Log all errors. */
    RXKB_LOG_LEVEL_WARNING = 30,  /**< Log warnings and errors. */
    RXKB_LOG_LEVEL_INFO = 40,     /**< Log information, warnings, and errors. */
    RXKB_LOG_LEVEL_DEBUG = 50     /**< Log everything. */
};

/**
 * Set the current logging level.
 *
 * @param ctx     The context in which to set the logging level.
 * @param level   The logging level to use.  Only messages from this level
 * and below will be logged.
 *
 * The default level is RXKB_LOG_LEVEL_ERROR.  The environment variable
 * RXKB_LOG_LEVEL, if set at the time the context was created, overrides the
 * default value.  It may be specified as a level number or name.
 */
void
rxkb_context_set_log_level(struct rxkb_context *ctx,
                           enum rxkb_log_level level);

/**
 * Get the current logging level.
 */
enum rxkb_log_level
rxkb_context_get_log_level(struct rxkb_context *ctx);

/**
 * Set a custom function to handle logging messages.
 *
 * @param ctx     The context in which to use the set logging function.
 * @param log_fn  The function that will be called for logging messages.
 * Passing NULL restores the default function, which logs to stderr.
 *
 * By default, log messages from this library are printed to stderr.  This
 * function allows you to replace the default behavior with a custom
 * handler.  The handler is only called with messages which match the
 * current logging level and verbosity settings for the context.
 * level is the logging level of the message.  @a format and @a args are
 * the same as in the vprintf(3) function.
 *
 * You may use rxkb_context_set_user_data() on the context, and then call
 * rxkb_context_get_user_data() from within the logging function to provide
 * it with additional private context.
 */
void
rxkb_context_set_log_fn(struct rxkb_context *ctx,
                        void (*log_fn)(struct rxkb_context *ctx,
                                       enum rxkb_log_level level,
                                       const char *format, va_list args));


/**
 * Parse the given ruleset. This can only be called once per context and once
 * parsed the data in the context is considered constant and will never
 * change.
 *
 * This function parses all files with the given ruleset name. See
 * rxkb_context_include_path_append() for details.
 *
 * If this function returns false, libxkbregistry failed to parse the xml files.
 * This is usually caused by invalid files on the host and should be debugged by
 * the host's administrator using external tools. Callers should reduce the
 * include paths to known good paths and/or fall back to a default RMLVO set.
 *
 * If this function returns false, the context should be be considered dead and
 * must be released with rxkb_context_unref().
 *
 * @param ctx The xkb registry context
 * @param ruleset The ruleset to parse, e.g. "evdev"
 * @return true on success or false on failure
 */
bool
rxkb_context_parse(struct rxkb_context *ctx, const char *ruleset);

/**
 * Parse the default ruleset as configured at build time. See
 * rxkb_context_parse() for details.
 */
bool
rxkb_context_parse_default_ruleset(struct rxkb_context *ctx);

/**
 * Increases the refcount of this object by one and returns the object.
 *
 * @param ctx The xkb registry context
 * @return The passed in object
 */
struct rxkb_context*
rxkb_context_ref(struct rxkb_context *ctx);

/**
 * Decreases the refcount of this object by one. Where the refcount of an
 * object hits zero, associated resources will be freed.
 *
 * @param ctx The xkb registry context
 * @return always NULL
 */
struct rxkb_context*
rxkb_context_unref(struct rxkb_context *ctx);

/**
 * Assign user-specific data. libxkbregistry will not look at or modify the
 * data, it will merely return the same pointer in
 * rxkb_context_get_user_data().
 *
 * @param ctx The xkb registry context
 * @param user_data User-specific data pointer
 */
void
rxkb_context_set_user_data(struct rxkb_context *ctx, void *user_data);

/**
 * Return the pointer passed into rxkb_context_get_user_data().
 *
 * @param ctx The xkb registry context
 * @return User-specific data pointer
 */
void *
rxkb_context_get_user_data(struct rxkb_context *ctx);

/**
 * Append a new entry to the context's include path.
 *
 * The include path handling is optimized for the most common use-case: a set of
 * system files that provide a complete set of MLVO and some
 * custom MLVO provided by a user **in addition** to the system set.
 *
 * The include paths should be given so that the least complete path is
 * specified first and the most complete path is appended last. For example:
 *
 * @code
 *    ctx = rxkb_context_new(RXKB_CONTEXT_NO_DEFAULT_INCLUDES);
 *    rxkb_context_include_path_append(ctx, "/home/user/.config/xkb");
 *    rxkb_context_include_path_append(ctx, "/usr/share/X11/xkb");
 *    rxkb_context_parse(ctx, "evdev");
 * @endcode
 *
 * The above example reflects the default behavior unless @ref
 * RXKB_CONTEXT_NO_DEFAULT_INCLUDES is provided.
 *
 * Loading of the files is in **reverse order**, i.e. the last path appended is
 * loaded first - in this case the ``/usr/share/X11/xkb`` path.
 * Any models, layouts, variants and options defined in the "evdev" ruleset
 * are loaded into the context. Then, any RMLVO found in the "evdev" ruleset of
 * the user's path (``/home/user/.config/xkb`` in this example) are **appended**
 * to the existing set.
 *
 * Note that data from previously loaded include paths is never overwritten,
 * only appended to. It is not not possible to change the system-provided data,
 * only to append new models, layouts, variants and options to it.
 *
 * In other words, to define a new variant of the "us" layout called "banana",
 * the following XML is sufficient.
 *
 * @verbatim
 * <xkbConfigRegistry version="1.1">
 * <layoutList>
 *   <layout>
 *     <configItem>
 *       <name>us</name>
 *     </configItem>
 *     <variantList>
 *       <variant>
 *         <configItem>
 *          <name>banana</name>
 *          <description>English (Banana)</description>
 *        </configItem>
 *      </variant>
 *    </layout>
 * </layoutList>
 * </xkbConfigRegistry>
 * @endverbatim
 *
 * The list of models, options and all other layouts (including "us" and its
 * variants) is taken from the system files. The resulting list of layouts will
 * thus have a "us" keyboard layout with the variant "banana" and all other
 * system-provided variants (dvorak, colemak, intl, etc.)
 *
 * This function must be called before rxkb_context_parse() or
 * rxkb_context_parse_default_ruleset().
 *
 * @returns true on success, or false if the include path could not be added
 * or is inaccessible.
 */
bool
rxkb_context_include_path_append(struct rxkb_context *ctx, const char *path);

/**
 * Append the default include paths to the context's include path.
 * See rxkb_context_include_path_append() for details about the merge order.
 *
 * This function must be called before rxkb_context_parse() or
 * rxkb_context_parse_default_ruleset().
 *
 * @returns true on success, or false if the include path could not be added
 * or is inaccessible.
 */
bool
rxkb_context_include_path_append_default(struct rxkb_context *ctx);

/**
 * Return the first model for this context. Use this to start iterating over
 * the models, followed by calls to rxkb_model_next(). Models are not sorted.
 *
 * The refcount of the returned model is not increased. Use rxkb_model_ref() if
 * you need to keep this struct outside the immediate scope.
 *
 * @return The first model in the model list.
 */
struct rxkb_model *
rxkb_model_first(struct rxkb_context *ctx);

/**
 * Return the next model for this context. Returns NULL when no more models
 * are available.
 *
 * The refcount of the returned model is not increased. Use rxkb_model_ref() if
 * you need to keep this struct outside the immediate scope.
 *
 * @return the next model or NULL at the end of the list
 */
struct rxkb_model *
rxkb_model_next(struct rxkb_model *m);

/**
 * Increase the refcount of the argument by one.
 *
 * @returns The argument passed in to this function.
 */
struct rxkb_model *
rxkb_model_ref(struct rxkb_model *m);

/**
 * Decrease the refcount of the argument by one. When the refcount hits zero,
 * all memory associated with this struct is freed.
 *
 * @returns always NULL
 */
struct rxkb_model *
rxkb_model_unref(struct rxkb_model *m);

/**
 * Return the name of this model. This is the value for M in RMLVO, to be used
 * with libxkbcommon.
 */
const char *
rxkb_model_get_name(struct rxkb_model *m);

/**
 * Return a human-readable description of this model. This function may return
 * NULL.
 */
const char *
rxkb_model_get_description(struct rxkb_model *m);

/**
 * Return the vendor name for this model. This function may return NULL.
 */
const char *
rxkb_model_get_vendor(struct rxkb_model *m);

/**
 * Return the popularity for this model.
 */
enum rxkb_popularity
rxkb_model_get_popularity(struct rxkb_model *m);

/**
 * Return the first layout for this context. Use this to start iterating over
 * the layouts, followed by calls to rxkb_layout_next(). Layouts are not sorted.
 *
 * The refcount of the returned layout is not increased. Use rxkb_layout_ref() if
 * you need to keep this struct outside the immediate scope.
 *
 * @return The first layout in the layout list.
 */
struct rxkb_layout *
rxkb_layout_first(struct rxkb_context *ctx);

/**
 * Return the next layout for this context. Returns NULL when no more layouts
 * are available.
 *
 * The refcount of the returned layout is not increased. Use rxkb_layout_ref()
 * if you need to keep this struct outside the immediate scope.
 *
 * @return the next layout or NULL at the end of the list
 */
struct rxkb_layout *
rxkb_layout_next(struct rxkb_layout *l);

/**
 * Increase the refcount of the argument by one.
 *
 * @returns The argument passed in to this function.
 */
struct rxkb_layout *
rxkb_layout_ref(struct rxkb_layout *l);

/**
 * Decrease the refcount of the argument by one. When the refcount hits zero,
 * all memory associated with this struct is freed.
 *
 * @returns always NULL
 */
struct rxkb_layout *
rxkb_layout_unref(struct rxkb_layout *l);

/**
 * Return the name of this layout. This is the value for L in RMLVO, to be used
 * with libxkbcommon.
 */
const char *
rxkb_layout_get_name(struct rxkb_layout *l);

/**
 * Return the variant of this layout. This is the value for V in RMLVO, to be
 * used with libxkbcommon.
 *
 * A variant does not stand on its own, it always depends on the base layout.
 * e.g. there may be multiple variants called "intl" but there is only one
 * "us(intl)".
 *
 * Where the variant is NULL, the layout is the base layout (e.g. "us").
 */
const char *
rxkb_layout_get_variant(struct rxkb_layout *l);

/**
 * Return a short (one-word) description of this layout. This function may
 * return NULL.
 */
const char *
rxkb_layout_get_brief(struct rxkb_layout *l);

/**
 * Return a human-readable description of this layout. This function may return
 * NULL.
 */
const char *
rxkb_layout_get_description(struct rxkb_layout *l);

/**
 * Return the popularity for this layout.
 */
enum rxkb_popularity
rxkb_layout_get_popularity(struct rxkb_layout *l);

/**
 * Return the first option group for this context. Use this to start iterating
 * over the option groups, followed by calls to rxkb_option_group_next().
 * Option groups are not sorted.
 *
 * The refcount of the returned option group is not increased. Use
 * rxkb_option_group_ref() if you need to keep this struct outside the immediate
 * scope.
 *
 * @return The first option group in the option group list.
 */
struct rxkb_option_group *
rxkb_option_group_first(struct rxkb_context *ctx);

/**
 * Return the next option group for this context. Returns NULL when no more
 * option groups are available.
 *
 * The refcount of the returned option group is not increased. Use
 * rxkb_option_group_ref() if you need to keep this struct outside the immediate
 * scope.
 *
 * @return the next option group or NULL at the end of the list
 */
struct rxkb_option_group *
rxkb_option_group_next(struct rxkb_option_group *g);

/**
 * Increase the refcount of the argument by one.
 *
 * @returns The argument passed in to this function.
 */
struct rxkb_option_group *
rxkb_option_group_ref(struct rxkb_option_group *g);

/**
 * Decrease the refcount of the argument by one. When the refcount hits zero,
 * all memory associated with this struct is freed.
 *
 * @returns always NULL
 */
struct rxkb_option_group *
rxkb_option_group_unref(struct rxkb_option_group *g);

/**
 * Return the name of this option group. This is **not** the value for O in
 * RMLVO, the name can be used for internal sorting in the caller. This function
 * may return NULL.
 */
const char *
rxkb_option_group_get_name(struct rxkb_option_group *m);

/**
 * Return a human-readable description of this option group. This function may
 * return NULL.
 */
const char *
rxkb_option_group_get_description(struct rxkb_option_group *m);

/**
 * @return true if multiple options within this option group can be selected
 *              simultaneously, false if all options within this option group
 *              are mutually exclusive.
 */
bool
rxkb_option_group_allows_multiple(struct rxkb_option_group *g);

/**
 * Return the popularity for this option group.
 */
enum rxkb_popularity
rxkb_option_group_get_popularity(struct rxkb_option_group *g);

/**
 * Return the first option for this option group. Use this to start iterating
 * over the options, followed by calls to rxkb_option_next(). Options are not
 * sorted.
 *
 * The refcount of the returned option is not increased. Use rxkb_option_ref()
 * if you need to keep this struct outside the immediate scope.
 *
 * @return The first option in the option list.
 */
struct rxkb_option *
rxkb_option_first(struct rxkb_option_group *group);

/**
 * Return the next option for this option group. Returns NULL when no more
 * options are available.
 *
 * The refcount of the returned options is not increased. Use rxkb_option_ref()
 * if you need to keep this struct outside the immediate scope.
 *
 * @returns The next option or NULL at the end of the list
 */
struct rxkb_option *
rxkb_option_next(struct rxkb_option *o);

/**
 * Increase the refcount of the argument by one.
 *
 * @returns The argument passed in to this function.
 */
struct rxkb_option *
rxkb_option_ref(struct rxkb_option *o);

/**
 * Decrease the refcount of the argument by one. When the refcount hits zero,
 * all memory associated with this struct is freed.
 *
 * @returns always NULL
 */
struct rxkb_option *
rxkb_option_unref(struct rxkb_option *o);

/**
 * Return the name of this option. This is the value for O in RMLVO, to be used
 * with libxkbcommon.
 */
const char *
rxkb_option_get_name(struct rxkb_option *o);

/**
 * Return a short (one-word) description of this option. This function may
 * return NULL.
 */
const char *
rxkb_option_get_brief(struct rxkb_option *o);

/**
 * Return a human-readable description of this option. This function may return
 * NULL.
 */
const char *
rxkb_option_get_description(struct rxkb_option *o);

/**
 * Return the popularity for this option.
 */
enum rxkb_popularity
rxkb_option_get_popularity(struct rxkb_option *o);

/**
 * Increase the refcount of the argument by one.
 *
 * @returns The argument passed in to this function.
 */
struct rxkb_iso639_code *
rxkb_iso639_code_ref(struct rxkb_iso639_code *iso639);

/**
 * Decrease the refcount of the argument by one. When the refcount hits zero,
 * all memory associated with this struct is freed.
 *
 * @returns always NULL
 */
struct rxkb_iso639_code *
rxkb_iso639_code_unref(struct rxkb_iso639_code *iso639);

/**
 * Return the ISO 639-3 code for this code (e.g. "eng", "fra").
 */
const char *
rxkb_iso639_code_get_code(struct rxkb_iso639_code *iso639);

/**
 * Return the first ISO 639 for this layout. Use this to start iterating over
 * the codes, followed by calls to rxkb_iso639_code_next(). Codes are not
 * sorted.
 *
 * The refcount of the returned code is not increased. Use rxkb_iso639_code_ref()
 * if you need to keep this struct outside the immediate scope.
 *
 * @return The first code in the code list.
 */
struct rxkb_iso639_code *
rxkb_layout_get_iso639_first(struct rxkb_layout *layout);

/**
 * Return the next code in the list. Returns NULL when no more codes
 * are available.
 *
 * The refcount of the returned codes is not increased. Use
 * rxkb_iso639_code_ref() if you need to keep this struct outside the immediate
 * scope.
 *
 * @returns The next code or NULL at the end of the list
 */
struct rxkb_iso639_code *
rxkb_iso639_code_next(struct rxkb_iso639_code *iso639);

/**
 * Increase the refcount of the argument by one.
 *
 * @returns The argument passed in to this function.
 */
struct rxkb_iso3166_code *
rxkb_iso3166_code_ref(struct rxkb_iso3166_code *iso3166);

/**
 * Decrease the refcount of the argument by one. When the refcount hits zero,
 * all memory associated with this struct is freed.
 *
 * @returns always NULL
 */
struct rxkb_iso3166_code *
rxkb_iso3166_code_unref(struct rxkb_iso3166_code *iso3166);

/**
 * Return the ISO 3166 Alpha 2 code for this code (e.g. "US", "FR").
 */
const char *
rxkb_iso3166_code_get_code(struct rxkb_iso3166_code *iso3166);

/**
 * Return the first ISO 3166 for this layout. Use this to start iterating over
 * the codes, followed by calls to rxkb_iso3166_code_next(). Codes are not
 * sorted.
 *
 * The refcount of the returned code is not increased. Use
 * rxkb_iso3166_code_ref() if you need to keep this struct outside the immediate
 * scope.
 *
 * @return The first code in the code list.
 */
struct rxkb_iso3166_code *
rxkb_layout_get_iso3166_first(struct rxkb_layout *layout);

/**
 * Return the next code in the list. Returns NULL when no more codes
 * are available.
 *
 * The refcount of the returned codes is not increased. Use
 * rxkb_iso3166_code_ref() if you need to keep this struct outside the immediate
 * scope.
 *
 * @returns The next code or NULL at the end of the list
 */
struct rxkb_iso3166_code *
rxkb_iso3166_code_next(struct rxkb_iso3166_code *iso3166);

/** @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _XKBREGISTRY_H_ */
