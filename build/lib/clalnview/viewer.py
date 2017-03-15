#!/usr/bin/env python

"""View commands."""

import curses
import click
import pickle
from .alphabet_colors import color_pairs
from .alterations import old_status, ref_map, pattern_ranges, trim_map
from os import path
from Bio import AlignIO
from collections import OrderedDict


# Some fixed variables
helpdict = OrderedDict([
    ('b', 'Phobicity based colors'),
    ('a', 'Normal view of amino acids'),
    ('z', 'Size based colors'),
    ('p', 'Previous file'),
    ('n', 'Next file'),
    ('o', 'Normal alignment'),
    ('l', 'Open location box'),
    ('f', 'Open pattern search box'),
    ('t', 'Alignment with mismaches only'),
    ('ctrl+up arrow', 'One page up'),
    ('ctrl+down arrow', 'One page down'),
    ('ctrl+left arrow', 'One page left'),
    ('click on left', 'Previous alignment'),
    ('click on right', 'Next alignment'),
    ('click on top', 'One line up'),
    ('click at bottom', 'One line down'),
    ('Click on sequence id', ' Reference based alignment'),
    ('ctrl+ rightarrow', 'One page right'),
    ('left arrow', '1 positions left'),
    ('right arrow', '1 positions right'),
    ('left arrow', '1 positions left'),
    ('up arrow', '1 positions up'),
    ('down arrow', '1 positions down'),
    ('End', 'End of alignment'),
    ('Home', 'Begining of alignment'),
    ('?', 'This help')
])


def is_alignment(alnfile):
    """Check for existance of alignemnt file."""
    for fltype in ['fasta', 'clustal', 'nexus', 'phylip']:
        try:
            alignment = AlignIO.read(open(alnfile), fltype)
            sequences = {}
            # alignment_length = alignment.get_alignment_length()
            for rec in alignment:
                sequences[rec.id] = rec.seq  # [:80]
            return sequences  # alignment_length
        except:
            continue
    else:
        return False


def show_me_the_alignment(file_list, hiddenpath, alntype):
    """Main display function."""
    screen = curses.initscr()
    curses.start_color()
    curses.use_default_colors()
    color_pair, alternative_color_pair = color_pairs(curses, "normal")
    current_x, current_y = 0, 0  # Left Top
    id_seq_gap = 15
    key_pressed = 100000

    current_file_num = 0
    try:
        current_file_num = int(open("%s/current_file_number" % hiddenpath).
                               read())
        if current_file_num > len(file_list):
            current_file_num = 0
    except:
        pass
    current_file = file_list[current_file_num]
    original = is_alignment(current_file)
    aln_file_found = False
    if original:
        aln_file_found = True
        sequence_ids = list(original.keys())  # If order changes,correct it
        sequence_ids.sort()
        indexes = None
        display_memory = old_status(current_file)
        current_x, current_y = (display_memory['display_x'],
                                display_memory['display_y'])
        if ((current_y > len(sequence_ids)) or
                (current_x > len(original[sequence_ids[0]]))):
            # To avoid buggy alteration of sequence files
            current_x, current_y = 0

        (color_pair,
         alternative_color_pair) = color_pairs(curses,
                                               display_memory["display_color"])
        if display_memory['display_mode'] == 't':
            key_pressed = ord('t')
        elif display_memory['display_mode'] == 'r':
            ref_name = display_memory['display_ref']
            sequences, sequence_ids = ref_map(original, ref_name)
            sequence_ids = sequence_ids[1:]
        else:
            sequences = original.copy()
            key_pressed = ord('o')
    else:
        key_pressed = ord('n')


    # Look up
    helpme = ['%15s: %s' % (k, helpdict[k]) for k in helpdict]
    helpme_height = len(helpme)
    helpme_width = max(map(len, helpme))

    curses.noecho()
    curses.curs_set(0)
    screen.keypad(1)
    curses.mousemask(1)
    curses.cbreak()
    max_y, max_x = screen.getmaxyx()
    vertical_shift = 0
    while key_pressed != ord('q'):
        # Handeling Resize of screen
        if key_pressed == curses.KEY_RESIZE:
            max_y, max_x = screen.getmaxyx()
            if max_y < 5 or max_x <= id_seq_gap + 1:
                screen.clear()
                key_pressed = screen.getch()
                continue
            if ((current_x + max_x - id_seq_gap - 1) >
                    len(sequences[sequence_ids[0]])):
                # Horizontal rearrangement
                if max_x - id_seq_gap - 1 < len(sequences[sequence_ids[0]]):
                    current_x = (len(sequences[sequence_ids[0]]) -
                                 (max_x - id_seq_gap - 1))
                else:
                    current_x = 0
            if ((current_y + max_y - vertical_shift - 2) >
                    len(sequence_ids)):
                current_y = (len(sequence_ids) -
                             (max_y - vertical_shift - 2))
                # Verical Rearrangement
                pass
        if alntype == "amn":
            if key_pressed == ord("b"):
                color_pair, alternative_color_pair = color_pairs(curses,
                                                                 "phobicity")
                display_memory["display_color"] = "phobicity"
            elif key_pressed == ord("z"):
                color_pair, alternative_color_pair = color_pairs(curses,
                                                                 "size")
                display_memory["display_color"] = "size"
            elif key_pressed == ord("a"):
                color_pair, alternative_color_pair = color_pairs(curses,
                                                                 "normal")
                display_memory["display_color"] = "normal"
        if key_pressed == ord('n') or key_pressed == ord('p'):

            if (key_pressed == ord('n') and
                    (current_file_num < len(file_list) - 1)):
                current_file_num += 1
            elif key_pressed == ord('p') and current_file_num > 0:
                current_file_num -= 1
            else:
                key_pressed = screen.getch()
                continue
            current_file = file_list[current_file_num]
            original = is_alignment(current_file)
            if not original:
                if (key_pressed == ord('n') and
                        (current_file_num < len(file_list) - 1)):
                    continue
                else:
                    if ((current_file_num == len(file_list) - 1) and
                            not aln_file_found):
                        break
                    else:
                        key_pressed = ord('p')
                        continue

                if current_file_num > 0:
                    # Less likely to be used unless some new files
                    # ToDo: Add a feature to remember the file number
                    if key_pressed == ord('p'):
                        continue
                else:
                    key_pressed == ord('n')

            sequence_ids = list(original.keys())
            sequences = original.copy()
            sequence_ids.sort()
            # If order changes,correct it
            display_memory = old_status(current_file)
            current_x, current_y = (display_memory['display_x'],
                                    display_memory['display_y'])
            if ((current_y > len(sequence_ids)) or
                    (current_x > len(original[sequence_ids[0]]))):
                # To avoid buggy alteration of sequence files
                current_x, current_y = 0
            (color_pair,
             alternative_color_pair) = color_pairs(curses,
                                                   display_memory[
                                                       "display_color"])

            if display_memory['display_mode'] == 't':
                # sequences, positions = trim_map(original)
                key_pressed = ord('t')
                continue

            elif display_memory['display_mode'] == 'r':

                ref_name = display_memory['display_ref']
                sequences, sequence_ids = ref_map(original, ref_name,
                                                  sequence_ids)
                sequence_ids = sequence_ids[1:]
                current_y -= vertical_shift
                vertical_shift = 1
            else:
                sequences = original.copy()
                key_pressed = ord('o')
                current_y -= vertical_shift
                vertical_shift = 1

        elif key_pressed == ord('l') and display_memory['display_mode'] != 't':
            search_box = screen.subwin(3, 30, int(max_y / 2), int(max_x / 2) -
                                       30)
            search_box.clear()
            search_box.box()
            in_str = "GoTo:"
            val_str = ""
            search_box.addstr(1, 1, in_str)
            search_box.refresh()
            curses.echo()
            char = search_box.getch()
            try:
                int(chr(char))
                val_str += chr(char)
            except:
                if char == 127:
                    if len(val_str) > 0:
                        val_str = val_str[:-1]

            while char != 27:  # Escape Key Value
                if char == 10:  # Enter Key Value
                    new_x = int(val_str)
                    break
                else:
                    search_box.clear()
                    search_box.box()
                    search_box.addstr(1, 1, in_str+val_str)
                    search_box.refresh()
                    char = search_box.getch()
                    try:
                        int(chr(char))
                        val_str += chr(char)
                    except:
                        if char == 127:
                            if len(val_str) > 0:
                                val_str = val_str[:-1]
            if char == 27:
                new_x = current_x + (max_x-id_seq_gap)/2

            search_box.clear()
            curses.noecho()

            if max_x - id_seq_gap >= len(sequences[sequence_ids[0]]):
                pass
            elif ((max_x - id_seq_gap)/2 <= new_x <=
                  (len(sequences[sequence_ids[0]]) - (max_x-id_seq_gap)/2)):
                current_x = int(new_x - (max_x - id_seq_gap)/2)
            elif (max_x - id_seq_gap)/2 > new_x:
                current_x = 0
            elif new_x > (len(sequences[sequence_ids[0]]) -
                          (max_x-id_seq_gap) / 2):
                current_x = int(len(sequences[sequence_ids[0]]) -
                                (max_x - id_seq_gap))

        elif (key_pressed == ord('f') and
              (display_memory['display_mode'] == 'o' or
               display_memory['display_mode'] == 'f')):
            ###
            search_box = screen.subwin(3, 50, int(max_y / 2), int(max_x / 2) -
                                       30)
            search_box.clear()
            search_box.box()
            in_str = "Pattern:"
            val_str = ""
            search_box.addstr(1, 1, in_str)
            search_box.refresh()
            curses.echo()
            char = search_box.getch()
            char_chr = chr(char)
            if char_chr.isalpha():
                val_str += char_chr.upper()

            while char != 27:  # Escape Key Value
                if char == 10:  # Enter Key Value
                    break
                else:
                    search_box.clear()
                    search_box.box()
                    search_box.addstr(1, 1, in_str+val_str)
                    search_box.refresh()
                    char = search_box.getch()
                    char_chr = chr(char)
                    if char_chr.isalpha():
                        val_str += char_chr.upper()
                    else:
                        if char == 127:
                            if len(val_str) > 0:
                                val_str = val_str[:-1]

            search_box.clear()
            curses.noecho()
            if char == 27 or len(val_str) == 0:
                key_pressed = 100000
                continue
            display_memory['display_mode'] = 'f'
            current_y -= vertical_shift
            vertical_shift = 0
            pattern_positions, search_movement = pattern_ranges(original,
                                                                val_str)
        #
        elif key_pressed == ord('o'):
            if display_memory['display_mode'] == 't':
                current_x = int(indexes[current_x].replace('.', ''))
            display_memory['display_mode'] = 'o'
            display_memory['display_ref'] = None
            sequences = original.copy()
            sequence_ids = list(sequences.keys())
            sequence_ids.sort()
            current_y -= vertical_shift
            vertical_shift = 0
        #
        elif key_pressed == ord('t'):
            display_memory['display_mode'] = 't'
            display_memory['display_ref'] = None
            sequences, indexes, max_height = trim_map(original)
            if current_x > len(indexes):
                current_x = len(indexes) - (max_x - id_seq_gap - 1)
                if current_x < 0:
                    current_x = 0
            current_y -= vertical_shift
            vertical_shift = max_height - 1

        elif key_pressed == ord('?'):
            xh, yh = 0, 0
            while key_pressed != ord('q'):
                if key_pressed == curses.KEY_RESIZE:
                    max_y, max_x = screen.getmaxyx()
                if key_pressed == curses.KEY_DOWN:
                    if yh + max_y == helpme_height-1:
                        key_pressed = screen.getch()
                        continue
                    yh += 1
                elif key_pressed == curses.KEY_UP:
                    if yh == 0:
                        key_pressed = screen.getch()
                        continue
                    yh -= 1

                elif key_pressed == curses.KEY_LEFT:
                    if xh == 0:
                        key_pressed = screen.getch()
                        continue
                    xh -= 1

                elif key_pressed == curses.KEY_RIGHT:
                    if xh + max_x == helpme_width or helpme_width <= max_x:
                        key_pressed = screen.getch()
                        continue
                    xh += 1
                screen.clear()
                for i, k in enumerate(helpme[yh: yh + max_y]):
                    screen.addstr(i, 0, k[xh: xh + max_x])
                key_pressed = screen.getch()
            key_pressed = 100000

        elif key_pressed == curses.KEY_MOUSE:

            _, mx, my, _, _ = curses.getmouse()
            if mx == 0:
                key_pressed = ord('p')
                continue
            elif mx == max_x - 1:
                key_pressed = ord('n')
                continue
            if my == 0:
                key_pressed = curses.KEY_UP
            elif my == max_y - 1:
                key_pressed = curses.KEY_DOWN

            if ((max_x < len(original[sequence_ids[0]]) +
                 id_seq_gap) and (mx > id_seq_gap)):
                if mx < (id_seq_gap + max_x) / 2:
                    if current_x < ((id_seq_gap + max_x) / 2 - mx):
                        shift = current_x
                    else:
                        shift = int((id_seq_gap + max_x) / 2 - mx)
                    current_x -= shift
                else:
                    if (len(original[sequence_ids[0]]) -
                        (current_x + max_x - id_seq_gap)
                        ) < (mx -
                             (id_seq_gap + max_x) / 2):
                        shift = int(len(original[sequence_ids[0]]) - (
                            current_x + max_x - id_seq_gap))
                    else:
                        shift = int(mx - (id_seq_gap + max_x) / 2)
                    current_x += shift

            if mx < 15 and my > 0:
                if display_memory['display_mode'] == 't':
                    if my < vertical_shift:
                        ref_name = sequence_ids[current_y]
                    else:
                        ref_name = sequence_ids[current_y + my -
                                                vertical_shift - 1]
                elif display_memory['display_mode'] == 'r':
                    if my < 2:
                        ref_name = sequence_ids[current_y]
                    else:
                        ref_name = sequence_ids[current_y + my - 2]
                else:
                    ref_name = sequence_ids[current_y + my - 1]
                sequences, sequence_ids = ref_map(original, ref_name)
                display_memory['display_mode'] = 'r'
                display_memory['display_ref'] = ref_name
                sequence_ids = sequence_ids[1:]
                current_y -= vertical_shift
                vertical_shift = 1

        elif key_pressed == 556:  # ctrl+right_key
            if current_x == (len(sequences[sequence_ids[0]]) -
                             (max_x - id_seq_gap - 1)):
                key_pressed = screen.getch()
                continue
            else:
                current_x += min([max_x - id_seq_gap,
                                  (len(sequences[sequence_ids[0]]) -
                                   (max_x - id_seq_gap +
                                    current_x - 1))])
        elif key_pressed == 541:  # ctrl+left_key
            if current_x == 0:
                key_pressed = screen.getch()
                continue
            else:
                if current_x <= max_x - id_seq_gap - 1:
                    current_x = 0
                else:
                    current_x -= max_x - id_seq_gap - 1
        elif key_pressed == 562:  # ctrl+up_key
            if current_y == 0:
                key_pressed = screen.getch()
                continue
            else:
                if current_y <= max_y - 2 - vertical_shift:
                    current_y = 0
                else:
                    current_y -= (max_y - 2 - vertical_shift)

        elif key_pressed == 521:  # ctrl+right_down
            if current_y == len(sequence_ids) - (max_y - 1 - vertical_shift):
                key_pressed = screen.getch()
                continue
            else:

                if (len(sequence_ids) - current_y) > 2*(max_y - 1 -
                                                        vertical_shift):
                    current_y += (max_y - 2 - vertical_shift)
                else:
                    current_y = len(sequence_ids) - (max_y - 2 -
                                                     vertical_shift)

        if key_pressed == curses.KEY_DOWN:
            if current_y == (len(sequence_ids) - (max_y - 2 - vertical_shift)):
                key_pressed = screen.getch()
                continue
            current_y += 1
        elif key_pressed == curses.KEY_UP:
            if current_y == 0:
                key_pressed = screen.getch()
                continue
            current_y -= 1
        elif key_pressed == curses.KEY_END:
            if current_x == (len(sequences[sequence_ids[0]]) -
                             (max_x - id_seq_gap - 1)):
                key_pressed = screen.getch()
                continue
            current_x = (len(sequences[sequence_ids[0]]) -
                         (max_x - id_seq_gap - 1))
        elif key_pressed == curses.KEY_HOME:
            if current_x == 0:
                key_pressed = screen.getch()
                continue
            current_x = 0
        elif key_pressed == curses.KEY_LEFT:
            if current_x == 0:
                key_pressed = screen.getch()
                continue
            else:
                current_x -= 1  # TODO:0 : Optimise it for error handling
        elif key_pressed == curses.KEY_RIGHT:
            if current_x == (len(sequences[sequence_ids[0]]) -
                             (max_x - id_seq_gap - 1)):
                key_pressed = screen.getch()
                continue
            elif ((current_x + max_x - id_seq_gap - 1) <
                  (len(sequences[sequence_ids[0]]))):
                current_x += 1

        screen.clear()
        if current_y < 0:
            current_y = 0

        if current_y + max_y - 2 - vertical_shift > len(sequence_ids):
            y_upto = len(sequence_ids)
        else:
            y_upto = current_y + max_y - 2 - vertical_shift
        if ((current_x + max_x - id_seq_gap - 1) >
                len(sequences[sequence_ids[0]])):
            x_upto = len(sequences[sequence_ids[0]])
        else:
            x_upto = current_x + max_x - id_seq_gap - 1

        if display_memory['display_mode'] == 'o':
            for i, s in enumerate(range(current_y, y_upto)):
                screen.addstr(i + 1, 1, sequence_ids[s])
                for j, n in enumerate(sequences[sequence_ids[s]
                                                ][current_x:x_upto]):
                    screen.addstr(i + 1, j + id_seq_gap, n, color_pair[n])
                    if i == 0 and (j + current_x) % 10 == 0:
                        if j + id_seq_gap > max_x - len(str(j+current_x)) - 1:
                            continue
                        screen.addstr(i, j + id_seq_gap,
                                      "|"+str(j + current_x))

        elif display_memory['display_mode'] == 'f':
            for i, s in enumerate(range(current_y, y_upto)):
                screen.addstr(i+1, 1, sequence_ids[s])
                for j, n in enumerate(sequences[sequence_ids[s]
                                                ][current_x:x_upto]):
                    if i == 0 and (j+current_x) % 10 == 0:
                            screen.addstr(i, j + id_seq_gap,
                                          "|"+str(j + current_x))
                    found = 0
                    for rng in pattern_positions[sequence_ids[s]]:
                        if rng[0] <= current_x + j < rng[1]:
                            found = 1
                            screen.addstr(i + 1, j + id_seq_gap, n,
                                          alternative_color_pair[n])
                    if not found:
                        screen.addstr(i + 1, j + id_seq_gap,
                                      n, color_pair[n])

        elif display_memory['display_mode'] == 't':
            for i, n in enumerate(indexes[current_x: x_upto]):

                for j, x in enumerate(n):
                    screen.addstr(j, i + id_seq_gap, x)
            for i, s in enumerate(range(current_y, y_upto),
                                  max_height):

                screen.addstr(i, 1, sequence_ids[s])
                for j, n in enumerate(sequences[sequence_ids[s]
                                                ][current_x: x_upto]):
                    screen.addstr(i, j + id_seq_gap, n, color_pair[n])

        elif display_memory['display_mode'] == 'r':
            screen.addstr(1, 1, ref_name)
            for j, n in enumerate(sequences[ref_name][current_x: x_upto]):
                if (j+current_x) % 10 == 0:
                    screen.addstr(0, j+id_seq_gap, "|"+str(j + current_x))
                screen.addstr(1, j+id_seq_gap, n, color_pair[n])
            for i, s in enumerate(range(current_y, y_upto)):
                screen.addstr(i + 2, 1, sequence_ids[s])

                for j, n in enumerate(sequences[sequence_ids[s]
                                                ][current_x: x_upto]):

                    screen.addstr(i + 2, j + id_seq_gap, n, color_pair[n])

        screen.addstr(max_y - 1, 0, current_file[:max_x - 1])
        screen.refresh()
        key_pressed = screen.getch()
        if chr(key_pressed) in 'npq':  # Next, Previous, Quit
            if display_memory['display_mode'] == 'f':
                display_memory['display_mode'] == 'o'
            display_memory['display_x'] = current_x
            display_memory['display_y'] = current_y
            hidden_path = path.split(current_file)
            hidden_path = "%s/.alnview/%s" % (hidden_path[0], hidden_path[1])
            pickle.dump(display_memory, open(hidden_path, "wb"))
            if aln_file_found:
                with open("%s/current_file_number" % hiddenpath,"w") as cfn:
                    cfn.write("%d" % current_file_num)

    curses.nocbreak()
    screen.keypad(0)
    curses.echo()
    curses.endwin()
    if not aln_file_found:
        click.echo("No valid alignment in supported"
                   " file format found. Exiting")