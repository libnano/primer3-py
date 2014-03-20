import os
import re
import warnings

from os.path import join as pjoin


def patchCfiles(package_dir, patch_fp):
    ''' Patch C files in a given package using a patch file with
    specially-delimited comments:

        Delimiter that prefaces an individual function / definition
            //#FILE#path/to/file/relative/to/source/package/dir##
            //#BLOCK#blockname##
            //#PRESIG#(?regex)*for[^pre]sig(usedtoidblockinoriginalfile)##
            //#POSTSIG#(?regex)*for[^post]sig(usedtoidblockinoriginalfile)##
            Content to insert
            //#ENDBLOCK##  

        Delimiter that prefaces a single statment to be replaced:
            //#FILE#path/to/file/relative/to/source/package/dir
            //#BLOCK#blockname
            //#REPLACE#(?regex)*for[^sig]tobe(replaced)
            //#STARTBLOCK##
            Content to insert
            //#ENDBLOCK##            

    The patched files are saved with the suffix _mod, so orig.c becomes
    orig_mod.c.
    '''

    def _processPatch(patched_content, fp):
        fp_split, ext = os.path.splitext(fp)
        fp = fp_split + '_mod' + ext
        with open(fp, 'w') as patch_fd:
            patch_fd.write(patched_content)

    def _correctHeaderRefs(patched_files):
        to_correct = {}
        for old_fp, new_fp in patched_files:
            directory, old_fn = os.path.split(old_fp)
            if old_fn.split('.')[-1] == 'h':
                _, new_fn = os.path.split(new_fp)
                to_correct.setdefault(directory, [])
                to_correct[directory].append((old_fn, new_fn))
        for directory, modded_fn_list in to_correct.items():
            files_to_check = os.listdir(directory)
            for fn in files_to_check:
                fn_pre, ext = os.path.splitext(fn)
                if fn_pre + '_mod' + ext in files_to_check:
                    continue
                if ext in ['.c', '.h', '.cpp']:
                    updated = False
                    with open(os.path.join(directory, fn)) as fd:
                        try:
                            contents = fd.read()
                            for old_fn, new_fn in modded_fn_list:
                                res = re.subn(
                                    '#include\s+[<|"]{}[>|"]'.format(old_fn),
                                    '#include "{}"'.format(new_fn),
                                    contents)
                                if res[1] > 0:
                                    updated = True
                                    contents = res[0]
                            if updated:
                                fn_base, ext = os.path.splitext(fn)
                                if not '_mod' in fn_base:
                                    fn_base += '_mod'
                                mod_fp = os.path.join(directory, fn_base + ext)
                                with open(mod_fp, 'w') as mod_fd:
                                    mod_fd.write(contents)
                                if ext == '.h':
                                    _correctHeaderRefs([(os.path.join(directory, fn), mod_fp)])
                        except UnicodeDecodeError:
                            pass

    _patch_re = re.compile(
        r'//#FILE#(?P<fp>.+?(?=##))##\s*'                                  \
        r'//#BLOCK#(?P<block>.+?(?=##))##\s*'                              \
        r'(?://#PRESIG#(?P<presig>.+?(?=##))##\s*)?'                       \
        r'(?://#POSTSIG#(?P<postsig>.+?(?=##))##\s*)?'                     \
        r'(?://#REPLACE#(?P<replace>.+?(?=##))##\s*)?'                     \
        r'//#STARTBLOCK##(?P<content>.+?(?=//#ENDBLOCK##))//#ENDBLOCK##',
        re.DOTALL)

    patch_dict = {}
    patched_files = []

    with open(patch_fp) as patch_fd:
        patches = re.finditer(_patch_re, patch_fd.read())

    for patch in patches:
        fp = pjoin(package_dir, patch.group('fp'))
        patch_dict.setdefault(fp, {})
        patch_dict[fp][patch.group('block')] = {
                'presig': patch.group('presig'),
                'postsig': patch.group('postsig'),
                'content': patch.group('content'),
                'replace': patch.group('replace')
            }
    for file_to_patch, patches in patch_dict.items():
        with open(pjoin(package_dir, file_to_patch)) as fd:
            contents_to_patch = fd.read()
        for block_name, patch in patches.items():
            if patch['replace'] and patch['replace'] != '':
                contents_to_patch = re.subn(patch['replace'], 
                                            patch['content'],
                                            contents_to_patch)[0]
            else:
                contents_to_patch = replaceCblock(contents_to_patch, 
                                                  patch['presig'],
                                                  patch['content'], 
                                                  patch['postsig'])
        _processPatch(contents_to_patch, file_to_patch)
        fp, ext = os.path.splitext(file_to_patch)
        new_fp = fp + '_mod' + ext
        patched_files.append((file_to_patch, new_fp))
    _correctHeaderRefs(patched_files)
    return patched_files


def replaceCblock(contents, presig, new_block, postsig=None):
    '''
    Replace a C block in string contents with `new_block`

    `presig` and `postsig` are a string, or preferably regular expressions, used
    to identify the bounds of the block to replace. replaceCblock will replace
    the entire block from the beginning of the `presig` to the end of `postsig`

    Raises ValueError if `presig` occurs multiple times or not at
    all and warns if `new_block` is already in place.

    ** TODO: implement postsig -- right now this just replaces the block
    starting with presig and ending when all brackets are balanced (also
    including the remainder of the last line with the final bracket)

    '''
    num_matches = len(re.findall(presig, contents))
    if num_matches > 1:
        raise ValueError('presig "{}" occurs {} times in contents'
                         .format(presig))
    elif contents.count(new_block):
        warnings.warn('New content already present in func with presig'
                      ' {}'.format(presig))
        return contents
    elif num_matches == 0:
        raise ValueError('Could not find presig "{}"'.format(presig))
    replace_start = ptr = re.search(presig, contents).start()
    while contents[ptr] != '{':
        ptr += 1
    left_brackets = 1
    right_brackets = 0
    while left_brackets != right_brackets:
        ptr += 1
        if contents[ptr] == '{':
            left_brackets += 1
        elif contents[ptr] == '}':
            right_brackets += 1
    while contents[ptr] not in ['\n', '\r']:
        ptr += 1
    replace_end = ptr + 1
    contents = contents[:replace_start] + new_block + contents[replace_end:]
    return contents
