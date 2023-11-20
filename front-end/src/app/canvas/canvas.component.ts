import { Component } from '@angular/core';
import { Block, BlockService } from '../block.service';

@Component({
  selector: 'app-canvas',
  templateUrl: './canvas.component.html',
  styleUrls: ['./canvas.component.css']
})

export class CanvasComponent {

  blockList: Block[] = []

  constructor(private blockService: BlockService) { 
    this.blockService.blockOnCanvas.subscribe(
      res => { this.blockList = res }
    )
  }

}
