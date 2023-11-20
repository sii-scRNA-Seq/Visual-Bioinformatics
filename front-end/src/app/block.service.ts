import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable, from } from 'rxjs';

@Injectable({
  providedIn: 'root'
})

export class BlockService {

  private readonly blockOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([])
  readonly blockOnCanvas: Observable<Block[]> = this.blockOnCanvas$.asObservable()

  constructor() { }

  addBlock(id: BlockId): void {
    switch(id) {

      case "LoadData": {
        if (this.blockOnCanvas$.getValue().length == 0) {
          this.blockOnCanvas$.next([{
            blockId: "LoadData",
            title: "Load Data", 
            possibleChildBlocks: [], 
            parameters: {},
            onRun: _ => from("")
          }])
        } else {
          console.log("You can't do that, it's wrong")
          /* this.blockOnCanvas$.next(
            this.blockOnCanvas$.getValue().concat([{
              title: "Load Data", 
              blockId: "LoadData", 
              onRun: _ => from(""), 
              possibleChildBlocks: [], 
              parameters: {}
            }])
          ) */
        }
        break;
      }


    }
  }

  removeBlock(id: BlockId): void {
    console.log(id);
    var newBlockList: Block[] = [];
    for (let i=0; i < this.blockOnCanvas$.getValue().length; i++) {
      if (this.blockOnCanvas$.getValue()[i].blockId == id) {
        this.blockOnCanvas$.next(newBlockList);
        break;
      } else {
        newBlockList.push(this.blockOnCanvas$.getValue()[i]);
      }
    }
    console.log("Block removed")
    // Remove the block that called this function

    // Either use blockId if each block can only exist once
    // Or add new position field and delete that block position

    // Could consider deleting only this block but that requires lots os logic
    // Easy solution is to delete this and all blocks below
  }

}

export interface Block {
  blockId: BlockId
  title: string
  possibleChildBlocks: BlockId[]
  parameters: Record<string, any>
  onRun: (block: Block) => Observable<any>
}

export type BlockId = "LoadData"