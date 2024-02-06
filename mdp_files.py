import glob
import os

if __name__ == "__main__":
    protocols = [os.path.basename(x).split('.')[0] for x in glob.glob('input/mdppath/*.X.mdp')]
    softcore = False

    for protocol in protocols:
        template = f'input/mdppath/{protocol}.X.mdp'

        with open(template, 'r') as f:
            content = f.read()

        for i in range(19):
            with open(f'input/mdppath/files/{protocol}.{i}.mdp', 'w') as f:
                if i <= 4:
                    content_new = content.replace('XXX', 'A')
                if i >= 5 and i <= 14:
                    content_new = content.replace('XXX', 'AB')
                elif i >= 15:
                    content_new = content.replace('XXX', 'B')
                if i ==0 or i==18:
                    content_new = content_new.replace('ZZZ', '500')
                else:
                    content_new = content_new.replace('ZZZ', '5000')
                content_new = content_new.replace('YYY', str(i))
                f.write(content_new)