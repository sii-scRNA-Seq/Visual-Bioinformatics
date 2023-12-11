import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable, from } from 'rxjs';
import { OutputService } from './output.service';
import { Parameter } from './parameter';

@Injectable({
  providedIn: 'root',
})

export class BlockService {
  private readonly blocksOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([]);
  readonly blocksOnCanvas: Observable<Block[]> = this.blocksOnCanvas$.asObservable();

  constructor(private outputService: OutputService) { }

  addBlock(id: BlockId): void {
    switch (id) {
      case 'loaddata': {
        if (this.blocksOnCanvas$.getValue().length == 0) {
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            title: 'Load Data',
            possibleChildBlocks: [],
            parameters: [],
            onRun: () => from(''),
          }]);
        }
        else {
          // Placeholder for action when block cannot be added
          console.log('You cant do that, its wrong');
        }
        break;
      }
      case 'basicfiltering': {  
        if (this.blocksOnCanvas$.getValue()[this.blocksOnCanvas$.getValue().length-1]?.blockId == 'loaddata') {
          const temp = this.blocksOnCanvas$.getValue();
          temp.push({
            blockId: 'basicfiltering',
            title: 'Basic Filtering',
            possibleChildBlocks: [],
            parameters: [
              {name: 'Minimum Cells Per Gene', value: 200},
              {name: 'Minimum Genes Per Cell', value: 3}
            ],
            onRun: () => from(''),
          });
          this.blocksOnCanvas$.next(temp);
        }
        else {
          // Placeholder for action when block cannot be added
          console.log('You cant do that, its wrong');
        }
        break;
      }
    }
  }

  removeBlock(id: BlockId): void {
    const newBlockList: Block[] = [];
    for (let i = 0; i < this.blocksOnCanvas$.getValue().length; i++) {
      if (this.blocksOnCanvas$.getValue()[i].blockId == id) {
        this.blocksOnCanvas$.next(newBlockList);
        break;
      }
      else {
        newBlockList.push(this.blocksOnCanvas$.getValue()[i]);
      }
    }
  }

  executeBlocks(): void {
    this.outputService.resetOutputs();
    for(let i=0; i < this.blocksOnCanvas$.getValue().length; i++) {
      this.outputService.executeBlock(this.blocksOnCanvas$.getValue()[i]);
    }
  }
}

export interface Block {
  blockId: BlockId
  title: string
  possibleChildBlocks: BlockId[]
  parameters: Parameter[]
  onRun: (block: Block) => Observable<unknown>
}

export type BlockId = 'loaddata' | 'basicfiltering';
