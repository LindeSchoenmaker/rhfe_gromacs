import glob
import os

if __name__ == "__main__":
    protocols = [os.path.basename(x).split('.')[0] for x in glob.glob('input/mdppath/*.X.mdp')]
    softcore = False

    for protocol in protocols:
        template = f'input/mdppath/{protocol}.X.mdp'

        with open(template, 'r') as f:
            content = f.read()

        if softcore:
            for i in range(11):
                with open(f'input/mdppath/files/{protocol}.{i}.mdp', 'w') as f:
                    content_new = content.replace('XXX', str(i))
                    f.write(content_new)

        else:
            for i in range(19):
                with open(f'input/mdppath/files/{protocol}.{i}.mdp', 'w') as f:
                    if i <= 4:
                        content_new = content.replace('XXX', 'A')
                    if i >= 5 and i <= 14:
                        content_new = content.replace('XXX', 'AB')
                    elif i >= 15:
                        content_new = content.replace('XXX', 'B')
                    content_new = content_new.replace('YYY', str(i))
                    f.write(content_new)