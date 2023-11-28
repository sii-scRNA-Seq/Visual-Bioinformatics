import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable, from } from 'rxjs';
import { OutputService } from './output.service';

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
            parameters: {},
            onRun: () => from(''),
          }]);
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
      const id = this.blocksOnCanvas$.getValue()[0].blockId;
      console.log(id);
      this.outputService.executeBlock(id);
    }
  }
}

export interface Block {
  blockId: BlockId
  title: string
  possibleChildBlocks: BlockId[]
  parameters: Record<string, unknown>
  onRun: (block: Block) => Observable<unknown>
}

export type BlockId = 'loaddata';
